/* Current calibrated dual-MTU model with external assistance.
 * At zero assistance, this is the verified BSE with standing feedback only.
 * Movement PD is gated by applied torque and corrects back activation.
 * Physical coefficients, Hill mechanics and angle-limit removal are unchanged. */
#include "../../../BSE/Direct_Activation_No_Limits/private/native/mtu_mechanics.h"
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#include <string.h>
#define NS 34
#define CS 69
#define NO 27

static double excitation(const double *unit,const double *p,int front,double q,double torque) {
    double lo=p[16+front]/unit[1],w=sin(p[23]);
    double slack=p[18+front]-lo*sqrt(unit[2]*unit[2]-w*w);
    double lmt=p[18+front]-(front?-1:1)*p[16+front]*(q-1.4);
    double ft=fmax(0,torque/unit[0]),lt=slack*(1+unit[6]*sqrt(ft));
    double x=(lmt-lt)/lo;
    if(x<=0)return 0;
    double l=sqrt(x*x+w*w),fl=exp(-pow((l-1)/.45,2));
    double a=(ft/(x/l)-passive(l))/fmax(fl,1e-8);
    return pow(clamp01(a),1/unit[5]);
}

static void standing_command(const double *p,const double *front,const double *cfg,
                             double q,double v,double mass,double integral,
                             double assist,double *command,double *request,double *correction) {
    double minimum[2];
    for(int u=0;u<2;u++) {
        double x,ft,ce,pe,vn,res;
        const double *unit=u?front:p;
        fiber(unit,q,0,0,.001,1,0,u?-1:1,p[16+u],p[18+u],sin(p[23]),&x,&ft,&ce,&pe,&vn,&res);
        minimum[u]=unit[0]*ft;
    }
    *correction=cfg[4]*tanh((cfg[1]*(cfg[0]-q)-cfg[2]*v+cfg[3]*integral)/cfg[4]);
    *request=-gravity(p,q,mass)+*correction-assist;
    double difference=*request-minimum[0]+minimum[1];
    command[0]=excitation(p,p,0,q,minimum[0]+cfg[5]+fmax(difference,0));
    command[1]=excitation(front,p,1,q,minimum[1]+cfg[5]+fmax(-difference,0));
}

static void output(const double *p,const double *s,double *o) {
    memset(o,0,NO*sizeof(double));
    o[0]=s[0];o[1]=s[1];o[2]=p[0]*s[15]-p[8]*s[16];
    o[3]=p[0]*s[15];o[4]=p[8]*s[16];
    o[5]=p[0]*s[17];o[6]=p[0]*s[19];o[7]=p[0]*.03*s[21]*s[23]/s[4];
    o[8]=p[8]*s[18];o[9]=p[8]*s[20];o[10]=p[8]*.03*s[22]*s[24]/s[5];
    o[11]=s[4];o[12]=s[5];o[13]=s[2];o[14]=s[3];
    o[15]=fmax(fabs(s[25]),fabs(s[26]));o[16]=s[27];
    o[17]=s[13];o[18]=s[14];o[19]=s[6];o[20]=s[11];o[21]=s[12];
    o[22]=s[7];o[26]=gravity(p,s[0],s[32]);
    /* The controlled wrapper fills target/assistance columns. Boundary is zero. */
}

static void initialize(const double *p,const double *cfg,const double *a,double mass,
                       double wait,const double *initial,double *s,double *o) {
    memset(s,0,NS*sizeof(double));
    double front[8]={p[8],p[10],p[11],0,p[20],p[9],p[21],p[22]},cmd[2]={a[0],a[1]};
    s[0]=initial[0];s[1]=initial[1];s[6]=wait>.5?1:0;
    if(s[6])standing_command(p,front,cfg,s[0],s[1],mass,0,0,cmd,&s[11],&s[12]);
    for(int u=0;u<2;u++) {
        const double *unit=u?front:p;
        s[13+u]=clamp01(cmd[u]);s[2+u]=pow(s[13+u],unit[5]);
        fiber(unit,s[0],s[2+u],0,.001,1,s[1],u?-1:1,p[16+u],p[18+u],sin(p[23]),
            &s[23+u],&s[15+u],&s[17+u],&s[19+u],&s[21+u],&s[25+u]);
        s[4+u]=sqrt(s[23+u]*s[23+u]+pow(sin(p[23]),2));
    }
    s[30]=a[0];s[31]=a[1];s[32]=mass;s[33]=wait;
    output(p,s,o);
}

static int advance_physical(const double *p,const double *cfg,double *s,const double *a,
                   double massEnd,double waitEnd,double assistEnd,const double *target,
                   const double *pd,double *filtered,double dt,int steps,double *o) {
    double front[8]={p[8],p[10],p[11],0,p[20],p[9],p[21],p[22]};
    double h=dt/steps,w=sin(p[23]);
    /* Release is instantaneous; an increase follows the external rise limit. */
    double startAssist=assistEnd<s[28]?assistEnd:s[28];
    for(int k=1;k<=steps;k++) {
        double f=(double)k/steps,mass=s[32]+f*(massEnd-s[32]);
        double assist=startAssist+f*(assistEnd-startAssist);
        double input[2]={clamp01(s[30]+f*(a[0]-s[30])),clamp01(s[31]+f*(a[1]-s[31]))};
        double cmd[2]={input[0],input[1]},request=0,correction=0;
        /* Standing feedback blends in and exits at the first action sample.
         * Movement activation PD is handled separately below. */
        if(waitEnd>.5) {
            s[6]=1+(s[6]-1)*exp(-h/cfg[6]);
            s[7]=fmax(-.5,fmin(.5,s[7]+h*(cfg[0]-s[0])));
            standing_command(p,front,cfg,s[0],s[1],mass,s[7],assist,cmd,&request,&correction);
        } else s[6]=0;
        /* Original movement PD: correct the BACK activation command only.
         * Targets are current causal unassisted-model states, never IK data.
         * The current measured assistance, not the requested command, gates PD. */
        if(waitEnd<=.5 && fabs(assist)>1e-9 && target && pd) {
            double delta=pd[0]*(target[2*(k-1)]-s[0])+pd[1]*(target[2*(k-1)+1]-s[1]);
            if(s[0]<1)delta=fmax(-pd[3],fmin(pd[2],delta));
            *filtered=pd[4]>0?delta+(*filtered-delta)*exp(-h/pd[4]):delta;
            input[0]=clamp01(input[0]+*filtered);
        } else *filtered=0;
        s[11]=s[6]*request;s[12]=s[6]*correction;
        for(int u=0;u<2;u++) {
            const double *unit=u?front:p;
            s[13+u]=s[6]*cmd[u]+(1-s[6])*input[u];
            double target=pow(clamp01(s[13+u]),unit[5]);
            s[2+u]=target+(s[2+u]-target)*exp(-h/unit[7]);
        }
        double inertia=mass*p[12]*p[12]+4*p[14]*p[15]*p[15]/3;
        double qold=s[0],vold=s[1],guess=qold+h*vold;
        double x[2],ft[2],ce[2],pe[2],vn[2],res[2];
        for(int iteration=0;iteration<9;iteration++) {
            for(int u=0;u<2;u++)fiber(u?front:p,guess,s[2+u],s[4+u],h,0,vold,u?-1:1,
                p[16+u],p[18+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
            double muscle=p[0]*ft[0]-p[8]*ft[1];
            double residual=inertia*(guess-qold-h*vold)/(h*h)+p[3]*(guess-qold)/h
                -muscle-gravity(p,guess,mass)-assist;
            if(fabs(residual)<1e-7)break;
            double xx,ff,cc,pp,vv,rr,eps=1e-6,plus=0;
            for(int u=0;u<2;u++) {
                fiber(u?front:p,guess+eps,s[2+u],s[4+u],h,0,vold,u?-1:1,
                    p[16+u],p[18+u],w,&xx,&ff,&cc,&pp,&vv,&rr);
                plus+=(u?-p[8]:p[0])*ff;
            }
            double derivative=inertia/(h*h)+p[3]/h-(plus-muscle)/eps
                -(gravity(p,guess+eps,mass)-gravity(p,guess,mass))/eps;
            /* This is a Newton-step bound, not a physical angle limit. */
            guess-=fmax(-.1,fmin(.1,residual/derivative));
        }
        s[0]=guess;s[1]=(guess-qold)/h;
        for(int u=0;u<2;u++) {
            fiber(u?front:p,guess,s[2+u],s[4+u],h,0,vold,u?-1:1,p[16+u],p[18+u],w,
                &s[23+u],&s[15+u],&s[17+u],&s[19+u],&s[21+u],&s[25+u]);
            s[4+u]=sqrt(s[23+u]*s[23+u]+w*w);
        }
        s[27]=inertia*(s[1]-vold)/h-p[0]*s[15]+p[8]*s[16]-gravity(p,s[0],mass)+p[3]*s[1]-assist;
        if(!isfinite(s[0]+s[1]+s[15]+s[16]) || fabs(s[25])+fabs(s[26])>1e-4 || fabs(s[27])>1)
            return 1;
    }
    s[28]=assistEnd;s[30]=a[0];s[31]=a[1];s[32]=massEnd;s[33]=waitEnd;
    output(p,s,o);return 0;
}


static void controlled_output(const double *p,const double *s,double *out) {
    output(p,s,out);out[24]=s[28];
    if(s[33]<=.5 && fabs(s[28])>1e-9)out[23]=s[34];
}

static void controlled_init(const double *p,const double *cfg,const double *a,double mass,
                            double wait,const double *initial,double *s,double *out) {
    memset(s,0,CS*sizeof(double));initialize(p,cfg,a,mass,wait,initial,s,out);
    memcpy(s+34,s,34*sizeof(double));controlled_output(p,s,out);
}

static int controlled_step(const double *p,const double *cfg,const double *pd,const double *actuator,
                           double *s,const double *a,double mass,double waiting,double command,
                           double dt,int steps,double *out) {
    double target[2048],row[NO],dummy=0;
    if(steps<1 || steps>1024)return 1;
    double applied=fmax(0,fmin(actuator[0],fmin(s[28]+actuator[1]*dt,command)));
    if(applied==0 && s[68]==0 && !memcmp(s,s+34,34*sizeof(double))) {
        /* Identical unassisted plant/shadow need only one physical solve. */
        if(advance_physical(p,cfg,s,a,mass,waiting,0,NULL,NULL,&dummy,dt,steps,row))return 1;
        memcpy(s+34,s,34*sizeof(double));controlled_output(p,s,out);return 0;
    }
    double start[2]={s[30],s[31]},startMass=s[32];
    /* Advance the causal zero-assistance shadow in lockstep. Each PD target
     * is the shadow state at the START of that same integration substep. */
    for(int j=0;j<steps;j++) {
        target[2*j]=s[34];target[2*j+1]=s[35];
        double f=(double)(j+1)/steps;
        double input[2]={start[0]+f*(a[0]-start[0]),start[1]+f*(a[1]-start[1])};
        if(advance_physical(p,cfg,s+34,input,startMass+f*(mass-startMass),waiting,
                            0,NULL,NULL,&dummy,dt/steps,1,row))return 1;
    }
    if(advance_physical(p,cfg,s,a,mass,waiting,applied,target,pd,s+68,dt,steps,row))return 1;
    controlled_output(p,s,out);return 0;
}

#ifdef MATLAB_MEX_FILE
static void require_vector(const mxArray *a,mwSize n) {
    if(!mxIsDouble(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=n)
        mexErrMsgIdAndTxt("CurrentMTU:Input","Unexpected real-double input size.");
    for(mwSize k=0;k<n;k++)if(!mxIsFinite(mxGetPr(a)[k]))
        mexErrMsgIdAndTxt("CurrentMTU:Input","Inputs must be finite.");
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    char action[16];
    if(nrhs<5 || mxGetString(prhs[0],action,sizeof(action)))
        mexErrMsgIdAndTxt("CurrentMTU:Call","Choose init, step or preview.");
    require_vector(prhs[1],24);require_vector(prhs[2],7);
    require_vector(prhs[3],5);require_vector(prhs[4],2);
    const double *p=mxGetPr(prhs[1]),*cfg=mxGetPr(prhs[2]),*pd=mxGetPr(prhs[3]),*actuator=mxGetPr(prhs[4]);
    if(!strcmp(action,"init")) {
        if(nrhs!=9 || nlhs!=2)mexErrMsgIdAndTxt("CurrentMTU:Call","Invalid initialization.");
        require_vector(prhs[5],2);require_vector(prhs[6],1);require_vector(prhs[7],1);require_vector(prhs[8],2);
        plhs[0]=mxCreateDoubleMatrix(CS,1,mxREAL);plhs[1]=mxCreateDoubleMatrix(1,NO,mxREAL);
        controlled_init(p,cfg,mxGetPr(prhs[5]),mxGetScalar(prhs[6]),mxGetScalar(prhs[7]),mxGetPr(prhs[8]),
                        mxGetPr(plhs[0]),mxGetPr(plhs[1]));
    } else if(!strcmp(action,"step")) {
        if(nrhs!=12 || nlhs!=2)mexErrMsgIdAndTxt("CurrentMTU:Call","Invalid advance.");
        require_vector(prhs[5],CS);require_vector(prhs[6],2);
        for(int k=7;k<12;k++)require_vector(prhs[k],1);
        double dt=mxGetScalar(prhs[10]);int steps=(int)mxGetScalar(prhs[11]);
        if(dt<=0 || steps<1)mexErrMsgIdAndTxt("CurrentMTU:Step","Positive step settings required.");
        plhs[0]=mxDuplicateArray(prhs[5]);plhs[1]=mxCreateDoubleMatrix(1,NO,mxREAL);
        if(controlled_step(p,cfg,pd,actuator,mxGetPr(plhs[0]),mxGetPr(prhs[6]),mxGetScalar(prhs[7]),
            mxGetScalar(prhs[8]),mxGetScalar(prhs[9]),dt,steps,mxGetPr(plhs[1])))
            mexErrMsgIdAndTxt("CurrentMTU:Equilibrium",
                "Hill equilibrium failed: q %.6g, v %.6g, momentum %.6g, fibers %.6g / %.6g.",
                mxGetPr(plhs[0])[0],mxGetPr(plhs[0])[1],mxGetPr(plhs[0])[27],
                mxGetPr(plhs[0])[25],mxGetPr(plhs[0])[26]);
    } else if(!strcmp(action,"preview")) {
        if(nrhs!=13 || nlhs!=1)mexErrMsgIdAndTxt("CurrentMTU:Call","Invalid preview.");
        require_vector(prhs[5],CS);mwSize n=mxGetM(prhs[7]),blocks=mxGetNumberOfElements(prhs[6]);
        if(n<1 || blocks<1 || mxGetN(prhs[7])!=2)mexErrMsgIdAndTxt("CurrentMTU:Input","Invalid forecast size.");
        require_vector(prhs[6],blocks);require_vector(prhs[7],2*n);require_vector(prhs[8],n);require_vector(prhs[9],n);
        require_vector(prhs[10],1);require_vector(prhs[11],1);require_vector(prhs[12],1);
        double dt=mxGetScalar(prhs[10]);int steps=(int)mxGetScalar(prhs[11]);
        double s[CS],row[NO];memcpy(s,mxGetPr(prhs[5]),sizeof(s));
        plhs[0]=mxCreateDoubleMatrix(n,NO,mxREAL);double *out=mxGetPr(plhs[0]);
        for(mwSize k=0;k<n;k++) {
            mwSize block=k*blocks/n;double input[2]={mxGetPr(prhs[7])[k],mxGetPr(prhs[7])[n+k]};
            if(controlled_step(p,cfg,pd,actuator,s,input,mxGetPr(prhs[8])[k],mxGetPr(prhs[9])[k],
                mxGetPr(prhs[6])[block],dt,steps,row)) {
                mxDestroyArray(plhs[0]);plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);return;
            }
            for(int j=0;j<NO;j++)out[k+n*j]=row[j];
        }
    } else mexErrMsgIdAndTxt("CurrentMTU:Call","Choose init, step or preview.");
}
#endif
