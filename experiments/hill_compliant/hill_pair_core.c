/* Compliant Hill muscle coupled to joint forward dynamics. No angle target.
 * Backward Euler is solved in the fiber and joint at every integration step.
 * Tendon force is monotone and tensile; no polynomial stiffness or force cap. */
#include <math.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
static double clamp01(double x){return fmax(0,fmin(1,x));}
static double fv(double v){return v<=0?fmax(0,1+v)/(1-v/.25):(1+1.5*v/.25)/(1+v/.25);}
static double passive(double l){return fmax(0,expm1(fmin(35,5*(l-1)/.6))/expm1(5));}
static void fiber(const double*p,double q,double a,double oldl,double h,int initial,double v0,int direction,
                  double R,double Lref,double w,
                  double *xn,double *ft,double *ce,double *pe,double *vn,double *res){
    (void)v0; /* Joint speed is independent of initial compliant-fiber speed. */
    double lo=R/p[1],xr=sqrt(p[2]*p[2]-w*w);
    double slack=Lref-lo*xr,c=p[6]*slack/lo,ell=xr-direction*p[1]*(q-1.4);
    /* Zero projected fiber length is the slack-tendon geometric limit.
     * A positive numerical floor would leave a spurious transmitted force. */
    /* A slack tendon may shorten below its slack length. The upper fiber
     * bound is the full available MT length, not only one tendon strain.
     * The previous ell+c bracket falsely excluded valid slack-tendon roots. */
    double lower=0,upper=fmax(1e-5,ell+slack/lo),x=initial?(lower+upper)/2:sqrt(fmax(0,oldl*oldl-w*w));
    x=fmax(lower+1e-10,fmin(upper-1e-10,x));
    for(int j=0;j<40;j++){
        double l=sqrt(x*x+w*w),cs=x/l,dcs=w*w/(l*l*l);
        /* A compliant tendon can change length while the initial fiber is
         * stationary. Initialize fiber equilibrium at zero fiber velocity;
         * the original joint velocity remains an independent initial state. */
        double vel=initial?0:(l-oldl)/(h*p[4]);
        double dvel=initial?0:cs/(h*p[4]);
        double stretch=fmax(0,ell-x),tendon=stretch*stretch/(c*c),dtendon=-2*stretch/(c*c);
        double fl=exp(-pow((l-1)/.45,2)),dfl=-2*(l-1)/(.45*.45)*fl*cs;
        double fvel=fv(vel),dfvel=vel<=-1?0:(vel<=0?5/pow(1-4*vel,2):2/pow(1+4*vel,2));
        double fp=passive(l),dfp=l<=1||5*(l-1)/.6>=35?0:(5/.6)*exp(5*(l-1)/.6)*cs/expm1(5);
        double active=a*fl*fvel,dactive=a*(dfl*fvel+fl*dfvel*dvel);
        double rhs=(active+fp+.03*vel)*cs;
        double f=tendon-rhs,der=dtendon-(dactive+dfp+.03*dvel)*cs-(active+fp+.03*vel)*dcs;
        if(fabs(f)<3e-11)break;
        if(f>0)lower=x;else upper=x;
        double next=x-f/der;
        if(!isfinite(next)||next<=lower||next>=upper)next=(lower+upper)/2;
        x=next;
    }
    double l=sqrt(x*x+w*w),cs=x/l;
    *vn=initial?0:(l-oldl)/(h*p[4]);
    *ft=pow(fmax(0,ell-x)/c,2); *ce=a*exp(-pow((l-1)/.45,2))*fv(*vn)*cs;
    *pe=passive(l)*cs; *xn=x; *res=*ft-*ce-*pe-.03*(*vn)*cs;
}
static double gravity(const double*p,double q,double m){return -9.81*p[18]*p[19]*cos(q)-9.81*m*p[12]*cos(q-p[13]);}
/* Elastic dead-zone soft limit. C1 onset, no angle clipping, no reference
 * feedback and no boundary damping. Its torque is NOT part of Hill torque. */
static double boundary(const double*p,double q){
    double low=fmax(p[30]-q,0),high=fmax(q-p[31],0);
    low=low*low/(low+.05);high=high*high/(high+.05);
    return p[29]*(tanh(p[28]*low/p[29])-tanh(p[28]*high/p[29]));
}
/* The two muscle units use identical Hill/tendon equations and opposing moment arms.
 * Output: q,v,netTorque,backTorque,frontTorque,backCE,backPE,backDamp,
 * frontCE,frontPE,frontDamp,backLength,frontLength,backActivation,frontActivation,
 * equilibriumResidual,momentumResidual. */
#if !defined(MATLAB_MEX_FILE) && defined(_WIN32)
__declspec(dllexport)
#elif !defined(MATLAB_MEX_FILE) && defined(__GNUC__)
__attribute__((visibility("default")))
#endif
int hill_pair_integrate(const double*p,const double*a,const double*m,int n,double q,double v,
                        double dt,int steps,double*out){
    double front[8]={p[8],p[10],p[11],0,p[24],p[9],p[25],p[26]};
    double h=dt/steps,w=sin(p[27]),x[2],ft[2],ce[2],pe[2],vn[2],res[2],oldl[2];
    double act[2]={pow(clamp01(a[0]),p[5]),pow(clamp01(a[n]),p[9])};
    for(int u=0;u<2;u++){
        const double*par=u==0?p:front;
        fiber(par,q,act[u],0,h,1,v,u==0?1:-1,p[20+u],p[22+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
        oldl[u]=sqrt(x[u]*x[u]+w*w);
    }
    int failed=0;
    for(int k=0;k<n;k++){
        double mom=0;
        if(k>0)for(int s=1;s<=steps;s++){
            double f=(double)s/steps,mass=m[k-1]+f*(m[k]-m[k-1]);
            for(int u=0;u<2;u++){
                const double*par=u==0?p:front;
                double input=a[u*n+k-1]+f*(a[u*n+k]-a[u*n+k-1]);
                double target=pow(clamp01(input),par[5]);
                act[u]=target+(act[u]-target)*exp(-h/par[7]);
            }
            double inertia=mass*p[12]*p[12]+4*p[18]*p[19]*p[19]/3;
            double qold=q,vold=v,guess=q+h*v;
            for(int j=0;j<9;j++){
                for(int u=0;u<2;u++)fiber(u==0?p:front,guess,act[u],oldl[u],h,0,vold,u==0?1:-1,p[20+u],p[22+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
                double muscle=p[0]*ft[0]-p[8]*ft[1];
                double fun=inertia*(guess-qold-h*vold)/(h*h)+p[3]*(guess-qold)/h-muscle-gravity(p,guess,mass)-boundary(p,guess);
                if(fabs(fun)<1e-7)break;
                double xx,ff,cc,pp,vv,rr,eps=1e-6,tplus=0;
                for(int u=0;u<2;u++){
                    fiber(u==0?p:front,guess+eps,act[u],oldl[u],h,0,vold,u==0?1:-1,p[20+u],p[22+u],w,&xx,&ff,&cc,&pp,&vv,&rr);
                    tplus+=(u==0?p[0]:-p[8])*ff;
                }
                double der=inertia/(h*h)+p[3]/h-(tplus-muscle)/eps-(gravity(p,guess+eps,mass)-gravity(p,guess,mass))/eps-(boundary(p,guess+eps)-boundary(p,guess))/eps;
                guess-=fmax(-.1,fmin(.1,fun/der));
            }
            q=guess;v=(q-qold)/h;
            for(int u=0;u<2;u++){
                fiber(u==0?p:front,q,act[u],oldl[u],h,0,vold,u==0?1:-1,p[20+u],p[22+u],w,&x[u],&ft[u],&ce[u],&pe[u],&vn[u],&res[u]);
                oldl[u]=sqrt(x[u]*x[u]+w*w);
            }
            mom=inertia*(v-vold)/h-p[0]*ft[0]+p[8]*ft[1]-gravity(p,q,mass)-boundary(p,q)+p[3]*v;
            if(!isfinite(q+v+ft[0]+ft[1])||fabs(q)>4||fabs(v)>50||fabs(res[0])+fabs(res[1])>1e-4||fabs(mom)>1){failed=1;break;}
        }
        out[k]=q;out[n+k]=v;out[2*n+k]=p[0]*ft[0]-p[8]*ft[1];
        out[3*n+k]=p[0]*ft[0];out[4*n+k]=p[8]*ft[1];
        out[5*n+k]=p[0]*ce[0];out[6*n+k]=p[0]*pe[0];out[7*n+k]=p[0]*.03*vn[0]*x[0]/oldl[0];
        out[8*n+k]=p[8]*ce[1];out[9*n+k]=p[8]*pe[1];out[10*n+k]=p[8]*.03*vn[1]*x[1]/oldl[1];
        out[11*n+k]=oldl[0];out[12*n+k]=oldl[1];out[13*n+k]=act[0];out[14*n+k]=act[1];
        out[15*n+k]=fmax(fabs(res[0]),fabs(res[1]));out[16*n+k]=mom;
        if(failed){for(int j=k+1;j<n;j++){out[j]=4;for(int col=1;col<17;col++)out[col*n+j]=0;}break;}
    }
    return failed;
}
#if defined(MATLAB_MEX_FILE) && !defined(HILL_CORE_NO_MEX_GATEWAY)
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    if(nrhs!=6||nlhs!=2)mexErrMsgIdAndTxt("hill:args","Six inputs and two outputs required.");
    for(int j=0;j<6;j++)if(!mxIsDouble(prhs[j])||mxIsComplex(prhs[j])||mxIsSparse(prhs[j]))mexErrMsgIdAndTxt("hill:type","Real full doubles required.");
    int n=(int)mxGetM(prhs[1]);
    if(mxGetNumberOfElements(prhs[0])!=32||mxGetN(prhs[1])!=2||mxGetNumberOfElements(prhs[2])!=(mwSize)n||mxGetNumberOfElements(prhs[3])!=2||n<2)mexErrMsgIdAndTxt("hill:size","Invalid input sizes (32 model parameters required).");
    plhs[0]=mxCreateDoubleMatrix(n,17,mxREAL); const double*y=mxGetPr(prhs[3]);
    int failed=hill_pair_integrate(mxGetPr(prhs[0]),mxGetPr(prhs[1]),mxGetPr(prhs[2]),n,y[0],y[1],mxGetScalar(prhs[4]),(int)mxGetScalar(prhs[5]),mxGetPr(plhs[0]));
    plhs[1]=mxCreateLogicalScalar(failed);
}
#endif
