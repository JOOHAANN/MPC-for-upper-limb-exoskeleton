/* Multisignal intent variant: the same Hill mechanics and neural excitation.
 * Additional fixed coefficients use abdominal drive and known load only.
 * EMG-informed motion intent and proprioception act only through excitation.
 * Angles and torques are outputs of the SAME compliant Hill forward system.
 * No reference trajectory, extra joint PD torque or per-condition parameter
 * switch enters this solver. This is a neural-control model, not open-loop
 * inference of motion from measured EMG alone. */
#define HILL_SEQUENCE_NO_MEX_GATEWAY
#include "hill_sequence_core.c"

/* The divisive candidate attenuates load-related tonic activation without
 * erasing a complete loaded movement. Both sources use the same decoder.
 * The default branch preserves the archived subtractive model exactly. */
static double decode_intent_signal(double input,double lead,double tonic,double gain){
#if defined(HILL_SOFT_TONIC_INTENT)
    double positive=fmax(0,input+lead),compensation=gain*tonic;
    if(compensation<=0)return clamp01(positive);
    return clamp01(positive*positive/(positive+compensation));
#elif defined(HILL_DIVISIVE_INTENT)
    return clamp01((input+lead)/(1+gain*tonic));
#else
    return clamp01(input+lead-gain*tonic);
#endif
}

/* Fixed cross-muscle/load features; no reference or trial identity. */
static double intent_input(double back,double front,double mass,double difference,const double*r){
    double load=mass/15.0;
    return back+r[10]*front+r[11]*front*load+r[12]*back*load
        +r[13]*back*back+r[14]*difference*load;
}

static void movement_command(const double*p,const double*front,const double*r,
                             double q,double v,double mass,double intent,
                             double front_input,double*command,double*desired,double*correction){
    double minimum[2],base_front;
    for(int u=0;u<2;u++){
        double x,ft,ce,pe,vn,res;
        const double*unit=u==0?p:front;
        fiber(unit,q,0,0,.001,1,0,u==0?1:-1,p[20+u],p[22+u],sin(p[27]),
              &x,&ft,&ce,&pe,&vn,&res);
        minimum[u]=unit[0]*ft;
    }
    double x,ft,ce,pe,vn,res;
    fiber(front,q,pow(clamp01(front_input),front[5]),0,.001,1,0,-1,p[21],p[23],sin(p[27]),
          &x,&ft,&ce,&pe,&vn,&res);
    base_front=front[0]*ft;
    *correction=r[3]*(intent-q)-r[4]*v;
    double request=-gravity(p,q,mass)+*correction;
    /* Asymmetric request saturation permits gravity-driven flexion with a
     * modest antagonist request. It does not cap actual muscle torque. */
    *desired=request<0?r[5]*tanh(request/r[5]):request;
    double back_target=fmax(minimum[0],*desired+base_front);
    double front_target=fmax(base_front,minimum[0]-*desired);
    back_target+=r[7];front_target+=r[7];
    command[0]=excitation_for_torque(p,p,0,q,back_target);
    command[1]=excitation_for_torque(front,p,1,q,front_target);
}

#if !defined(MATLAB_MEX_FILE) && defined(_WIN32)
__declspec(dllexport)
#elif !defined(MATLAB_MEX_FILE) && defined(__GNUC__)
__attribute__((visibility("default")))
#endif
int hill_neural_integrate(const double*p,const double*a,const double*m,
                         const double*wait,int n,const double*initial,
                         double dt,int steps,const double*cfg,const double*r,double*out){
    double front[8]={p[8],p[10],p[11],0,p[24],p[9],p[25],p[26]};
    double q=initial[0],v=initial[1],h=dt/steps,w=sin(p[27]);
    double x[2],ft[2],ce[2],pe[2],vn[2],res[2],oldl[2],act[2],command[2];
    double blend=clamp01(wait[0]),integral=0,desired=0,correction=0,control[2]={0,0};
    double fast=a[0],slow=a[0];
    /* Only the intent decoder removes the estimated load-related tonic
     * increment. This subtraction does not modify the source input itself.
     * A zero coefficient exactly reproduces the earlier decoder. */
    double neutral=excitation_for_torque(p,p,0,cfg[0],fmax(0,-gravity(p,cfg[0],0)));
    double tonic=0;
    if(r[9]>0)tonic=fmax(0,excitation_for_torque(p,p,0,cfg[0],fmax(0,-gravity(p,cfg[0],m[0])))-neutral);
    double intent=cfg[0]+r[15]-r[0]*pow(decode_intent_signal(intent_input(a[0],a[n],m[0],0,r),0,tonic,r[9]),r[1]);
    if(blend>1e-8)stance_command(p,front,cfg,q,v,m[0],integral,control,&desired,&correction);
    else movement_command(p,front,r,q,v,m[0],intent,a[n],control,&desired,&correction);
    for(int u=0;u<2;u++){
        const double*unit=u==0?p:front;
        double feedback=blend+(1-blend)*r[6];
        command[u]=(1-feedback)*clamp01(a[u*n])+feedback*control[u];
        act[u]=pow(clamp01(command[u]),unit[5]);
        fiber(unit,q,act[u],0,h,1,v,u==0?1:-1,p[20+u],p[22+u],w,
              &x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
        oldl[u]=sqrt(x[u]*x[u]+w*w);
    }
    int failed=0;
    for(int k=0;k<n;k++){
        double mom=0;
        if(k>0)for(int s=1;s<=steps;s++){
            double fraction=(double)s/steps;
            double mass=m[k-1]+fraction*(m[k]-m[k-1]);
            double target_blend=clamp01(wait[k-1]+fraction*(wait[k]-wait[k-1]));
            blend=target_blend+(blend-target_blend)*exp(-h/cfg[6]);
            integral=fmax(-.5,fmin(.5,integral+h*blend*(cfg[0]-q)));
            double input[2];
            for(int u=0;u<2;u++)input[u]=clamp01(a[u*n+k-1]+fraction*(a[u*n+k]-a[u*n+k-1]));
            /* Causal lead from EMG history compensates filtering and motor
             * delay; no future sample or reference kinematics is used. */
            fast=input[0]+(fast-input[0])*exp(-h/.06);
            slow=input[0]+(slow-input[0])*exp(-h/.25);
            tonic=0;
            if(r[9]>0)tonic=fmax(0,excitation_for_torque(p,p,0,cfg[0],fmax(0,-gravity(p,cfg[0],mass)))-neutral);
            double intent_signal=decode_intent_signal(intent_input(input[0],input[1],mass,fast-slow,r),r[8]*(fast-slow)/(.25-.06),tonic,r[9]);
            double target_intent=cfg[0]+r[15]-r[0]*pow(intent_signal,r[1]);
            intent=target_intent+(intent-target_intent)*exp(-h/r[2]);
            double posture[2]={0,0},moving[2]={0,0},dpost=0,cpost=0,dmove=0,cmove=0;
            if(blend>1e-8)stance_command(p,front,cfg,q,v,mass,integral,posture,&dpost,&cpost);
            if(blend<1-1e-8)movement_command(p,front,r,q,v,mass,intent,input[1],moving,&dmove,&cmove);
            desired=blend*dpost+(1-blend)*dmove;
            correction=blend*cpost+(1-blend)*cmove;
            for(int u=0;u<2;u++){
                const double*unit=u==0?p:front;
                command[u]=blend*posture[u]+(1-blend)*((1-r[6])*input[u]+r[6]*moving[u]);
                double target=pow(clamp01(command[u]),unit[5]);
                act[u]=target+(act[u]-target)*exp(-h/unit[7]);
            }
            /* Identical fiber/tendon/joint solve to hill_sequence_integrate. */
            double inertia=mass*p[12]*p[12]+4*p[18]*p[19]*p[19]/3;
            double qold=q,vold=v,guess=q+h*v;
            for(int j=0;j<9;j++){
                for(int u=0;u<2;u++)fiber(u==0?p:front,guess,act[u],oldl[u],h,0,vold,u==0?1:-1,
                    p[20+u],p[22+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
                double muscle=p[0]*ft[0]-p[8]*ft[1];
                double fun=inertia*(guess-qold-h*vold)/(h*h)+p[3]*(guess-qold)/h-muscle-gravity(p,guess,mass)-boundary(p,guess);
                if(fabs(fun)<1e-7)break;
                double xx,ff,cc,pp,vv,rr,eps=1e-6,tplus=0;
                for(int u=0;u<2;u++){
                    fiber(u==0?p:front,guess+eps,act[u],oldl[u],h,0,vold,u==0?1:-1,
                          p[20+u],p[22+u],w,&xx,&ff,&cc,&pp,&vv,&rr);
                    tplus+=(u==0?p[0]:-p[8])*ff;
                }
                double der=inertia/(h*h)+p[3]/h-(tplus-muscle)/eps
                           -(gravity(p,guess+eps,mass)-gravity(p,guess,mass))/eps
                           -(boundary(p,guess+eps)-boundary(p,guess))/eps;
                guess-=fmax(-.1,fmin(.1,fun/der));
            }
            q=guess;v=(q-qold)/h;
            for(int u=0;u<2;u++){
                fiber(u==0?p:front,q,act[u],oldl[u],h,0,vold,u==0?1:-1,
                      p[20+u],p[22+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
                oldl[u]=sqrt(x[u]*x[u]+w*w);
            }
            mom=inertia*(v-vold)/h-p[0]*ft[0]+p[8]*ft[1]-gravity(p,q,mass)-boundary(p,q)+p[3]*v;
            if(!isfinite(q+v+ft[0]+ft[1])||fabs(q)>4||fabs(v)>50
               ||fabs(res[0])+fabs(res[1])>1e-4||fabs(mom)>1){failed=1;break;}
        }
        out[k]=q;out[n+k]=v;out[2*n+k]=p[0]*ft[0]-p[8]*ft[1];
        out[3*n+k]=p[0]*ft[0];out[4*n+k]=p[8]*ft[1];
        out[5*n+k]=p[0]*ce[0];out[6*n+k]=p[0]*pe[0];out[7*n+k]=p[0]*.03*vn[0]*x[0]/oldl[0];
        out[8*n+k]=p[8]*ce[1];out[9*n+k]=p[8]*pe[1];out[10*n+k]=p[8]*.03*vn[1]*x[1]/oldl[1];
        out[11*n+k]=oldl[0];out[12*n+k]=oldl[1];out[13*n+k]=act[0];out[14*n+k]=act[1];
        out[15*n+k]=fmax(fabs(res[0]),fabs(res[1]));out[16*n+k]=mom;
        out[17*n+k]=command[0];out[18*n+k]=command[1];out[19*n+k]=blend;
        out[20*n+k]=desired;out[21*n+k]=correction;out[22*n+k]=integral;
        out[23*n+k]=intent;
        if(failed){for(int j=k+1;j<n;j++){out[j]=4;for(int col=1;col<24;col++)out[col*n+j]=0;}break;}
    }
    return failed;
}

#if defined(MATLAB_MEX_FILE) && !defined(HILL_NEURAL_NO_MEX_GATEWAY)
void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[]){
    if(nrhs!=9||nlhs!=2)mexErrMsgIdAndTxt("hillNeural:args","Nine inputs and two outputs required.");
    for(int j=0;j<9;j++)if(!mxIsDouble(prhs[j])||mxIsComplex(prhs[j])||mxIsSparse(prhs[j]))
        mexErrMsgIdAndTxt("hillNeural:type","Real full doubles required.");
    int n=(int)mxGetM(prhs[1]);
    if(mxGetNumberOfElements(prhs[0])!=32||mxGetN(prhs[1])!=2||mxGetNumberOfElements(prhs[2])!=(mwSize)n
       ||mxGetNumberOfElements(prhs[3])!=(mwSize)n||mxGetNumberOfElements(prhs[4])!=2
       ||mxGetNumberOfElements(prhs[7])!=7||mxGetNumberOfElements(prhs[8])!=16||n<2)
        mexErrMsgIdAndTxt("hillNeural:size","Invalid input sizes.");
    plhs[0]=mxCreateDoubleMatrix(n,24,mxREAL);
    int failed=hill_neural_integrate(mxGetPr(prhs[0]),mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),n,
        mxGetPr(prhs[4]),mxGetScalar(prhs[5]),(int)mxGetScalar(prhs[6]),mxGetPr(prhs[7]),mxGetPr(prhs[8]),mxGetPr(plhs[0]));
    plhs[1]=mxCreateLogicalScalar(failed);
}
#endif
