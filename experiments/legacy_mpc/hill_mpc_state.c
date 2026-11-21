/* Stateful access to the SELECTED neural/Hill model for predictive control.
 * The original fiber/tendon laws, calibrated parameters and human motion
 * decoder are reused. Exoskeleton torque enters joint momentum balance and
 * the human muscle request. No measured angle is an input to this plant. */
#define HILL_NEURAL_NO_MEX_GATEWAY
#include "../hill_compliant/hill_neural_core.c"
#include <string.h>
#if defined(_WIN32)
#define API __declspec(dllexport)
#else
#define API __attribute__((visibility("default")))
#endif
#define STATE_SIZE 34
#define OUTPUT_SIZE 27

static void assisted_request(const double*p,const double*front,const double*r,
                             const double*cfg,double q,double v,double mass,
                             double intent,double front_input,double integral,
                             double assist,int standing,double*cmd,double*desired,double*corr){
    if(assist==0){
        if(standing)stance_command(p,front,cfg,q,v,mass,integral,cmd,desired,corr);
        else movement_command(p,front,r,q,v,mass,intent,front_input,cmd,desired,corr);
        return;
    }
    double minimum[2];
    for(int u=0;u<2;u++){
        double x,ft,ce,pe,vn,res;
        const double*unit=u==0?p:front;
        fiber(unit,q,0,0,.001,1,0,u==0?1:-1,p[20+u],p[22+u],sin(p[27]),&x,&ft,&ce,&pe,&vn,&res);
        minimum[u]=unit[0]*ft;
    }
    if(standing){
        *corr=cfg[4]*tanh((cfg[1]*(cfg[0]-q)-cfg[2]*v+cfg[3]*integral)/cfg[4]);
        *desired=-gravity(p,q,mass)+*corr-assist;
        double difference=*desired-minimum[0]+minimum[1];
        cmd[0]=excitation_for_torque(p,p,0,q,minimum[0]+cfg[5]+fmax(difference,0));
        cmd[1]=excitation_for_torque(front,p,1,q,minimum[1]+cfg[5]+fmax(-difference,0));
    }else{
        double x,ft,ce,pe,vn,res;
        fiber(front,q,pow(clamp01(front_input),front[5]),0,.001,1,0,-1,p[21],p[23],sin(p[27]),&x,&ft,&ce,&pe,&vn,&res);
        double base_front=front[0]*ft;
        *corr=r[3]*(intent-q)-r[4]*v;
        double request=-gravity(p,q,mass)+*corr-assist;
        *desired=request<0?r[5]*tanh(request/r[5]):request;
        cmd[0]=excitation_for_torque(p,p,0,q,fmax(minimum[0],*desired+base_front)+r[7]);
        cmd[1]=excitation_for_torque(front,p,1,q,fmax(base_front,minimum[0]-*desired)+r[7]);
    }
}

static void state_output(const double*p,const double*s,double*out){
    out[0]=s[0];out[1]=s[1];out[2]=p[0]*s[15]-p[8]*s[16];
    out[3]=p[0]*s[15];out[4]=p[8]*s[16];
    out[5]=p[0]*s[17];out[6]=p[0]*s[19];out[7]=p[0]*.03*s[21]*s[23]/s[4];
    out[8]=p[8]*s[18];out[9]=p[8]*s[20];out[10]=p[8]*.03*s[22]*s[24]/s[5];
    out[11]=s[4];out[12]=s[5];out[13]=s[2];out[14]=s[3];
    out[15]=fmax(fabs(s[25]),fabs(s[26]));out[16]=s[27];
    out[17]=s[13];out[18]=s[14];out[19]=s[6];out[20]=s[11];out[21]=s[12];out[22]=s[7];
    out[23]=s[10];out[24]=s[28];out[25]=boundary(p,s[0]);out[26]=gravity(p,s[0],s[32]);
}

API int hill_mpc_init(const double*p,const double*r,const double*cfg,const double*a,
                      double mass,double waiting,const double*initial,double*s,double*out){
    memset(s,0,STATE_SIZE*sizeof(double));
    double front[8]={p[8],p[10],p[11],0,p[24],p[9],p[25],p[26]},cmd[2],desired,corr;
    s[0]=initial[0];s[1]=initial[1];s[6]=clamp01(waiting);s[8]=a[0];s[9]=a[0];
    s[29]=excitation_for_torque(p,p,0,cfg[0],fmax(0,-gravity(p,cfg[0],0)));
    double tonic=r[9]>0?fmax(0,excitation_for_torque(p,p,0,cfg[0],fmax(0,-gravity(p,cfg[0],mass)))-s[29]):0;
    s[10]=cfg[0]-r[0]*pow(decode_intent_signal(a[0],0,tonic,r[9]),r[1]);
    if(s[6]>1e-8)stance_command(p,front,cfg,s[0],s[1],mass,0,cmd,&desired,&corr);
    else movement_command(p,front,r,s[0],s[1],mass,s[10],a[1],cmd,&desired,&corr);
    s[11]=desired;s[12]=corr;
    for(int u=0;u<2;u++){
        const double*unit=u==0?p:front;
        double feedback=s[6]+(1-s[6])*r[6];
        s[13+u]=(1-feedback)*clamp01(a[u])+feedback*cmd[u];
        s[2+u]=pow(clamp01(s[13+u]),unit[5]);
        fiber(unit,s[0],s[2+u],0,.001,1,s[1],u==0?1:-1,p[20+u],p[22+u],sin(p[27]),
              &s[23+u],&s[15+u],&s[17+u],&s[19+u],&s[21+u],&s[25+u]);
        s[4+u]=sqrt(s[23+u]*s[23+u]+pow(sin(p[27]),2));
    }
    s[30]=a[0];s[31]=a[1];s[32]=mass;s[33]=waiting;
    state_output(p,s,out);
    return 0;
}

API int hill_mpc_advance(const double*p,const double*r,const double*cfg,double*s,
                         const double*a,double massEnd,double waitEnd,double assistEnd,
                         double dt,int steps,double*out){
    double front[8]={p[8],p[10],p[11],0,p[24],p[9],p[25],p[26]};
    double h=dt/steps,w=sin(p[27]),startU=s[28];
    for(int step=1;step<=steps;step++){
        double f=(double)step/steps,mass=s[32]+f*(massEnd-s[32]),assist=startU+f*(assistEnd-startU);
        double btarget=clamp01(s[33]+f*(waitEnd-s[33]));
        s[6]=btarget+(s[6]-btarget)*exp(-h/cfg[6]);
        s[7]=fmax(-.5,fmin(.5,s[7]+h*s[6]*(cfg[0]-s[0])));
        double input[2]={clamp01(s[30]+f*(a[0]-s[30])),clamp01(s[31]+f*(a[1]-s[31]))};
        s[8]=input[0]+(s[8]-input[0])*exp(-h/.06);
        s[9]=input[0]+(s[9]-input[0])*exp(-h/.25);
        double tonic=r[9]>0?fmax(0,excitation_for_torque(p,p,0,cfg[0],fmax(0,-gravity(p,cfg[0],mass)))-s[29]):0;
        double signal=decode_intent_signal(input[0],r[8]*(s[8]-s[9])/(.25-.06),tonic,r[9]);
        double targetIntent=cfg[0]-r[0]*pow(signal,r[1]);
        s[10]=targetIntent+(s[10]-targetIntent)*exp(-h/r[2]);
        double post[2]={0,0},moving[2]={0,0},dp=0,cp=0,dm=0,cm=0;
        if(s[6]>1e-8)assisted_request(p,front,r,cfg,s[0],s[1],mass,s[10],input[1],s[7],assist,1,post,&dp,&cp);
        if(s[6]<1-1e-8)assisted_request(p,front,r,cfg,s[0],s[1],mass,s[10],input[1],s[7],assist,0,moving,&dm,&cm);
        s[11]=s[6]*dp+(1-s[6])*dm;s[12]=s[6]*cp+(1-s[6])*cm;
        for(int u=0;u<2;u++){
            const double*unit=u==0?p:front;
            s[13+u]=s[6]*post[u]+(1-s[6])*((1-r[6])*input[u]+r[6]*moving[u]);
            double target=pow(clamp01(s[13+u]),unit[5]);
            s[2+u]=target+(s[2+u]-target)*exp(-h/unit[7]);
        }
        double inertia=mass*p[12]*p[12]+4*p[18]*p[19]*p[19]/3,qold=s[0],vold=s[1],guess=qold+h*vold;
        double x[2],ft[2],ce[2],pe[2],vn[2],res[2];
        for(int j=0;j<9;j++){
            for(int u=0;u<2;u++)fiber(u==0?p:front,guess,s[2+u],s[4+u],h,0,vold,u==0?1:-1,
                p[20+u],p[22+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
            double muscle=p[0]*ft[0]-p[8]*ft[1];
            double fun=inertia*(guess-qold-h*vold)/(h*h)+p[3]*(guess-qold)/h-muscle-gravity(p,guess,mass)-boundary(p,guess)-assist;
            if(fabs(fun)<1e-7)break;
            double xx,ff,cc,pp,vv,rr,eps=1e-6,tplus=0;
            for(int u=0;u<2;u++){
                fiber(u==0?p:front,guess+eps,s[2+u],s[4+u],h,0,vold,u==0?1:-1,
                      p[20+u],p[22+u],w,&xx,&ff,&cc,&pp,&vv,&rr);
                tplus+=(u==0?p[0]:-p[8])*ff;
            }
            double der=inertia/(h*h)+p[3]/h-(tplus-muscle)/eps-(gravity(p,guess+eps,mass)-gravity(p,guess,mass))/eps
                       -(boundary(p,guess+eps)-boundary(p,guess))/eps;
            guess-=fmax(-.1,fmin(.1,fun/der));
        }
        s[0]=guess;s[1]=(guess-qold)/h;
        for(int u=0;u<2;u++){
            fiber(u==0?p:front,guess,s[2+u],s[4+u],h,0,vold,u==0?1:-1,p[20+u],p[22+u],w,
                  &s[23+u],&s[15+u],&s[17+u],&s[19+u],&s[21+u],&s[25+u]);
            s[4+u]=sqrt(s[23+u]*s[23+u]+w*w);
        }
        s[27]=inertia*(s[1]-vold)/h-p[0]*s[15]+p[8]*s[16]-gravity(p,s[0],mass)-boundary(p,s[0])+p[3]*s[1]-assist;
        if(!isfinite(s[0]+s[1]+s[15]+s[16])||fabs(s[0])>4||fabs(s[1])>50
           ||fabs(s[25])+fabs(s[26])>1e-4||fabs(s[27])>1)return 1;
    }
    s[28]=assistEnd;s[30]=a[0];s[31]=a[1];s[32]=massEnd;s[33]=waitEnd;
    state_output(p,s,out);
    return 0;
}

/* Causal hold of current EMG, load and wait state; no future measurements. */
API int hill_mpc_predict(const double*p,const double*r,const double*cfg,const double*initial,
                         const double*blocks,int blockCount,int points,double dt,int steps,
                         double rate,double*out){
    double s[STATE_SIZE];memcpy(s,initial,sizeof(s));
    double input[2]={s[30],s[31]},mass=s[32],waiting=s[33];
    for(int k=0;k<points;k++){
        int block=k*blockCount/points;
        double target=blocks[block],u=fmax(s[28]-rate*dt,fmin(s[28]+rate*dt,target));
        if(hill_mpc_advance(p,r,cfg,s,input,mass,waiting,u,dt,steps,out+k*OUTPUT_SIZE))return 1;
    }
    return 0;
}

/* Causal, bounded input-trend forecast. The slope comes only from past EMG.
 * Extrapolation decays over tau seconds instead of growing over the horizon. */
API int hill_mpc_predict_input_trend(const double*p,const double*r,const double*cfg,
                         const double*initial,const double*blocks,const double*slope,
                         int blockCount,int points,double dt,int steps,double rate,
                         double tau,double*out){
    double s[STATE_SIZE];memcpy(s,initial,sizeof(s));
    double input0[2]={s[30],s[31]},mass=s[32],waiting=s[33];
    for(int k=0;k<points;k++){
        int block=k*blockCount/points;
        double target=blocks[block],u=fmax(s[28]-rate*dt,fmin(s[28]+rate*dt,target));
        double elapsed=(k+1)*dt,gain=tau*(1-exp(-elapsed/tau));
        double input[2]={clamp01(input0[0]+slope[0]*gain),clamp01(input0[1]+slope[1]*gain)};
        if(hill_mpc_advance(p,r,cfg,s,input,mass,waiting,u,dt,steps,out+k*OUTPUT_SIZE))return 1;
    }
    return 0;
}

/* Forecast from a trained activation generator and a declared task schedule.
 * The caller must not populate these arrays with future measured EMG or IK.
 * This entry point changes prediction only; the actual Hill plant is unchanged. */
API int hill_mpc_predict_map(const double*p,const double*r,const double*cfg,
                         const double*initial,const double*blocks,
                         const double*activation,const double*mass,const double*waiting,
                         int blockCount,int points,double dt,int steps,double rate,double*out){
    double s[STATE_SIZE];memcpy(s,initial,sizeof(s));
    for(int k=0;k<points;k++){
        int block=k*blockCount/points;
        double target=blocks[block],u=fmax(s[28]-rate*dt,fmin(s[28]+rate*dt,target));
        if(hill_mpc_advance(p,r,cfg,s,activation+2*k,mass[k],waiting[k],u,dt,steps,out+k*OUTPUT_SIZE))return 1;
    }
    return 0;
}
