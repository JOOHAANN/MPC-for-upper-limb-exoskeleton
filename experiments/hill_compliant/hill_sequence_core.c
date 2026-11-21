/* Continuous three-load rollout with a muscle-driven standing controller.
 * The fiber/tendon/joint equations are shared with hill_pair_core.c.
 * The stance PD acts through excitation, never as an extra balance torque.
 * One initialization per sequence; no state reset at a load or wait boundary. */
#define HILL_CORE_NO_MEX_GATEWAY
#include "hill_pair_core.c"

/* Static force inversion is used ONLY to request neural commands.
 * Actual forces always come from the compliant Hill solve below. */
static double excitation_for_torque(const double *unit, const double *p,
                                   int index, double q, double torque) {
    double lo=p[20+index]/unit[1],w=sin(p[27]);
    double slack=p[22+index]-lo*sqrt(unit[2]*unit[2]-w*w);
    double lmt=p[22+index]-(index==0?1:-1)*p[20+index]*(q-1.4);
    double ft=fmax(0,torque/unit[0]);
    double lt=slack*(1+unit[6]*sqrt(ft));
    double x=(lmt-lt)/lo;
    if(x<=0)return 0; /* No feasible tensile force at this slack geometry. */
    double l=sqrt(x*x+w*w),cs=x/l,fl=exp(-pow((l-1)/.45,2));
    double activation=(ft/cs-passive(l))/fmax(fl,1e-8);
    return pow(clamp01(activation),1/unit[5]);
}

static void stance_command(const double *p, const double *front, const double *cfg,
                           double q, double v, double mass, double integral,
                           double *command, double *desired, double *correction) {
    double passive_torque[2];
    for(int u=0;u<2;u++){
        double x,ft,ce,pe,vn,res;
        const double *unit=u==0?p:front;
        fiber(unit,q,0,0,.001,1,0,u==0?1:-1,p[20+u],p[22+u],sin(p[27]),
              &x,&ft,&ce,&pe,&vn,&res);
        passive_torque[u]=unit[0]*ft;
    }
    *correction=cfg[4]*tanh((cfg[1]*(cfg[0]-q)-cfg[2]*v+cfg[3]*integral)/cfg[4]);
    *desired=-gravity(p,q,mass)+*correction;
    double difference=*desired-passive_torque[0]+passive_torque[1];
    double target[2]={passive_torque[0]+cfg[5]+fmax(difference,0),
                      passive_torque[1]+cfg[5]+fmax(-difference,0)};
    command[0]=excitation_for_torque(p,p,0,q,target[0]);
    command[1]=excitation_for_torque(front,p,1,q,target[1]);
}

#if !defined(MATLAB_MEX_FILE) && defined(_WIN32)
__declspec(dllexport)
#elif !defined(MATLAB_MEX_FILE) && defined(__GNUC__)
__attribute__((visibility("default")))
#endif
int hill_sequence_integrate(const double*p,const double*a,const double*m,
                            const double*wait,int n,const double*initial,
                            double dt,int steps,const double*cfg,double*out){
    double front[8]={p[8],p[10],p[11],0,p[24],p[9],p[25],p[26]};
    double q=initial[0],v=initial[1],h=dt/steps,w=sin(p[27]);
    double x[2],ft[2],ce[2],pe[2],vn[2],res[2],oldl[2],act[2],command[2];
    double blend=clamp01(wait[0]),integral=0,desired=0,correction=0,control[2]={0,0};
    if(blend>1e-8)stance_command(p,front,cfg,q,v,m[0],integral,control,&desired,&correction);
    for(int u=0;u<2;u++){
        const double*unit=u==0?p:front;
        command[u]=(1-blend)*clamp01(a[u*n])+blend*control[u];
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
            if(blend<1e-8){blend=0;integral*=exp(-h);desired=0;correction=0;}
            if(blend>1e-8)stance_command(p,front,cfg,q,v,mass,integral,control,&desired,&correction);
            for(int u=0;u<2;u++){
                const double*unit=u==0?p:front;
                double input=a[u*n+k-1]+fraction*(a[u*n+k]-a[u*n+k-1]);
                command[u]=(1-blend)*clamp01(input)+blend*control[u];
                double target=pow(clamp01(command[u]),unit[5]);
                act[u]=target+(act[u]-target)*exp(-h/unit[7]);
            }
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
        if(failed){for(int j=k+1;j<n;j++){out[j]=4;for(int col=1;col<23;col++)out[col*n+j]=0;}break;}
    }
    return failed;
}

#if defined(MATLAB_MEX_FILE) && !defined(HILL_SEQUENCE_NO_MEX_GATEWAY)
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    if(nrhs!=8||nlhs!=2)mexErrMsgIdAndTxt("hillSequence:args","Eight inputs and two outputs required.");
    for(int j=0;j<8;j++)if(!mxIsDouble(prhs[j])||mxIsComplex(prhs[j])||mxIsSparse(prhs[j]))
        mexErrMsgIdAndTxt("hillSequence:type","Real full doubles required.");
    int n=(int)mxGetM(prhs[1]);
    if(mxGetNumberOfElements(prhs[0])!=32||mxGetN(prhs[1])!=2||
       mxGetNumberOfElements(prhs[2])!=(mwSize)n||mxGetNumberOfElements(prhs[3])!=(mwSize)n||
       mxGetNumberOfElements(prhs[4])!=2||mxGetNumberOfElements(prhs[7])!=7||n<2)
        mexErrMsgIdAndTxt("hillSequence:size","Invalid input sizes.");
    plhs[0]=mxCreateDoubleMatrix(n,23,mxREAL);
    int failed=hill_sequence_integrate(mxGetPr(prhs[0]),mxGetPr(prhs[1]),mxGetPr(prhs[2]),
        mxGetPr(prhs[3]),n,mxGetPr(prhs[4]),mxGetScalar(prhs[5]),(int)mxGetScalar(prhs[6]),
        mxGetPr(prhs[7]),mxGetPr(plhs[0]));
    plhs[1]=mxCreateLogicalScalar(failed);
}
#endif
