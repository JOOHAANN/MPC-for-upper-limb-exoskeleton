/* Direct activation -> two compliant MTUs -> one joint.
 * Standing balance is the only human feedback. No movement intent, movement
 * PD, angle boundary torque, angle clipping or reference trajectory exists. */
#include "mtu_mechanics.h"
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#include <string.h>
#define NS 34
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
                             double *command,double *request,double *correction) {
    double minimum[2];
    for(int u=0;u<2;u++) {
        double x,ft,ce,pe,vn,res;
        const double *unit=u?front:p;
        fiber(unit,q,0,0,.001,1,0,u?-1:1,p[16+u],p[18+u],sin(p[23]),&x,&ft,&ce,&pe,&vn,&res);
        minimum[u]=unit[0]*ft;
    }
    *correction=cfg[4]*tanh((cfg[1]*(cfg[0]-q)-cfg[2]*v+cfg[3]*integral)/cfg[4]);
    *request=-gravity(p,q,mass)+*correction;
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
    /* Columns 24/25/26 remain zero: no intent, assistance or angle boundary. */
}

static void initialize(const double *p,const double *cfg,const double *a,double mass,
                       double wait,const double *initial,double *s,double *o) {
    memset(s,0,NS*sizeof(double));
    double front[8]={p[8],p[10],p[11],0,p[20],p[9],p[21],p[22]},cmd[2]={a[0],a[1]};
    s[0]=initial[0];s[1]=initial[1];s[6]=wait>.5?1:0;
    if(s[6])standing_command(p,front,cfg,s[0],s[1],mass,0,cmd,&s[11],&s[12]);
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

static int advance(const double *p,const double *cfg,double *s,const double *a,
                   double massEnd,double waitEnd,double dt,int steps,double *o) {
    double front[8]={p[8],p[10],p[11],0,p[20],p[9],p[21],p[22]};
    double h=dt/steps,w=sin(p[23]);
    for(int k=1;k<=steps;k++) {
        double f=(double)k/steps,mass=s[32]+f*(massEnd-s[32]);
        double input[2]={clamp01(s[30]+f*(a[0]-s[30])),clamp01(s[31]+f*(a[1]-s[31]))};
        double cmd[2]={input[0],input[1]},request=0,correction=0;
        /* Standing can blend in smoothly. It exits immediately at an action
         * sample: zero posture feedback is guaranteed throughout movement. */
        if(waitEnd>.5) {
            s[6]=1+(s[6]-1)*exp(-h/cfg[6]);
            s[7]=fmax(-.5,fmin(.5,s[7]+h*(cfg[0]-s[0])));
            standing_command(p,front,cfg,s[0],s[1],mass,s[7],cmd,&request,&correction);
        } else s[6]=0;
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
                -muscle-gravity(p,guess,mass);
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
        s[27]=inertia*(s[1]-vold)/h-p[0]*s[15]+p[8]*s[16]-gravity(p,s[0],mass)+p[3]*s[1];
        if(!isfinite(s[0]+s[1]+s[15]+s[16]) || fabs(s[25])+fabs(s[26])>1e-4 || fabs(s[27])>1)
            return 1;
    }
    s[30]=a[0];s[31]=a[1];s[32]=massEnd;s[33]=waitEnd;
    output(p,s,o);return 0;
}

#if defined(_WIN32)
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif
/* Training and MATLAB replay call the same integration routines above. */
EXPORT int direct_mtu_rollout(const double *p,const double *cfg,int count,
                             const double *input,const double *mass,const double *waiting,
                             const double *initial,double dt,int steps,double *result) {
    double state[NS],row[NO];
    initialize(p,cfg,input,mass[0],waiting[0],initial,state,row);
    memcpy(result,row,sizeof(row));
    for(int k=1;k<count;k++) {
        if(advance(p,cfg,state,input+2*k,mass[k],waiting[k],dt,steps,row))return k;
        memcpy(result+NO*k,row,sizeof(row));
    }
    return count;
}

#ifdef MATLAB_MEX_FILE
static void vector(const mxArray *a,mwSize n) {
    if(!mxIsDouble(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=n)
        mexErrMsgIdAndTxt("DirectMTU:Input","Unexpected real-double input size.");
    for(mwSize j=0;j<n;j++)if(!mxIsFinite(mxGetPr(a)[j]))
        mexErrMsgIdAndTxt("DirectMTU:Input","Inputs must be finite.");
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    char action[16];
    if(nrhs<3 || mxGetString(prhs[0],action,sizeof(action)))
        mexErrMsgIdAndTxt("DirectMTU:Call","Choose init or step.");
    vector(prhs[1],24);vector(prhs[2],7);
    const double *p=mxGetPr(prhs[1]),*cfg=mxGetPr(prhs[2]);
    if(!strcmp(action,"rollout")) {
        if(nrhs!=9 || nlhs!=2)mexErrMsgIdAndTxt("DirectMTU:Call","Invalid rollout.");
        mwSize n=mxGetM(prhs[3]);
        if(n<1 || mxGetN(prhs[3])!=2)mexErrMsgIdAndTxt("DirectMTU:Input","Expected N-by-2 input.");
        vector(prhs[3],2*n);vector(prhs[4],n);vector(prhs[5],n);vector(prhs[6],2);
        vector(prhs[7],1);vector(prhs[8],1);
        double dt=mxGetScalar(prhs[7]);int steps=(int)mxGetScalar(prhs[8]);
        if(dt<=0 || steps<1)mexErrMsgIdAndTxt("DirectMTU:Step","Positive integration settings required.");
        double *inputs=mxMalloc(2*n*sizeof(double)),*rows=mxMalloc(NO*n*sizeof(double));
        for(mwSize k=0;k<n;k++){inputs[2*k]=mxGetPr(prhs[3])[k];inputs[2*k+1]=mxGetPr(prhs[3])[n+k];}
        for(mwSize k=0;k<NO*n;k++)rows[k]=mxGetNaN();
        int completed=direct_mtu_rollout(p,cfg,(int)n,inputs,mxGetPr(prhs[4]),mxGetPr(prhs[5]),
            mxGetPr(prhs[6]),dt,steps,rows);
        plhs[0]=mxCreateDoubleMatrix(n,NO,mxREAL);plhs[1]=mxCreateDoubleScalar(completed);
        for(mwSize k=0;k<n;k++)for(int j=0;j<NO;j++)mxGetPr(plhs[0])[k+n*j]=rows[NO*k+j];
        mxFree(inputs);mxFree(rows);
    } else if(!strcmp(action,"init")) {
        if(nrhs!=7 || nlhs!=2)mexErrMsgIdAndTxt("DirectMTU:Call","Invalid initialization.");
        vector(prhs[3],2);vector(prhs[4],1);vector(prhs[5],1);vector(prhs[6],2);
        plhs[0]=mxCreateDoubleMatrix(NS,1,mxREAL);plhs[1]=mxCreateDoubleMatrix(1,NO,mxREAL);
        initialize(p,cfg,mxGetPr(prhs[3]),mxGetScalar(prhs[4]),mxGetScalar(prhs[5]),
            mxGetPr(prhs[6]),mxGetPr(plhs[0]),mxGetPr(plhs[1]));
    } else if(!strcmp(action,"step")) {
        if(nrhs!=9 || nlhs!=2)mexErrMsgIdAndTxt("DirectMTU:Call","Invalid step.");
        vector(prhs[3],NS);vector(prhs[4],2);
        for(int j=5;j<9;j++)vector(prhs[j],1);
        double dt=mxGetScalar(prhs[7]);int steps=(int)mxGetScalar(prhs[8]);
        if(dt<=0 || steps<1)mexErrMsgIdAndTxt("DirectMTU:Step","Positive integration settings required.");
        plhs[0]=mxDuplicateArray(prhs[3]);plhs[1]=mxCreateDoubleMatrix(1,NO,mxREAL);
        if(advance(p,cfg,mxGetPr(plhs[0]),mxGetPr(prhs[4]),mxGetScalar(prhs[5]),
            mxGetScalar(prhs[6]),dt,steps,mxGetPr(plhs[1])))
            mexErrMsgIdAndTxt("DirectMTU:Equilibrium","Fiber/joint equilibrium failed; no output was clipped.");
    } else mexErrMsgIdAndTxt("DirectMTU:Call","Choose init or step.");
}
#endif
