/* Optional accelerator for hill_rollout.m. Equations and RK4 are identical.
 * Compile in MATLAB with: mex -R2017b hill_integrate_mex.c
 * The readable MATLAB implementation remains the fallback and test oracle. */
#include "mex.h"
#include <math.h>
static double acceleration(const double *p,double q,double v,double a,double m) {
    double w=p[1]*sin(.24086), z=p[1]*sqrt(p[3]*p[3]-pow(sin(.24086),2))-p[2]*(q-1.4);
    double l=sqrt(z*z+w*w), c=fmax(z,0)/l, vn=-p[2]*v*c/(p[5]*p[1]);
    double fl=exp(-pow((l/p[1]-1)/.45,2));
    double fv=vn<=0?fmax(0,1+vn)/(1-vn/.25):(1+1.5*vn/.25)/(1+vn/.25);
    double fp=fmax(0,(exp(fmin(40,5*(l/p[1]-1)/.6))-1)/(exp(5)-1));
    double torque=p[2]*p[0]*(pow(fmax(0,fmin(1,a)),p[6])*fl*fv+fp)*c;
    double grav=-9.81*22.6006*.345*cos(q)-9.81*m*.760606*cos(q-.936847);
    return (torque+grav-p[4]*v)/(m*.760606*.760606+4*22.6006*.345*.345/3);
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    if(nrhs!=6 || nlhs!=2) mexErrMsgIdAndTxt("hill:args","Six inputs and two outputs required.");
    for(int j=0;j<6;j++) if(!mxIsDouble(prhs[j])||mxIsComplex(prhs[j])||mxIsSparse(prhs[j]))
        mexErrMsgIdAndTxt("hill:type","Inputs must be real full doubles.");
    mwSize n=mxGetNumberOfElements(prhs[1]);
    if(mxGetNumberOfElements(prhs[0])!=7 || n<2 || mxGetNumberOfElements(prhs[2])!=n ||
       (mxGetNumberOfElements(prhs[3])!=1 && mxGetNumberOfElements(prhs[3])!=2) || mxGetNumberOfElements(prhs[4])!=1 || mxGetNumberOfElements(prhs[5])!=1)
        mexErrMsgIdAndTxt("hill:size","Invalid input sizes.");
    const double *p=mxGetPr(prhs[0]),*a=mxGetPr(prhs[1]),*m=mxGetPr(prhs[2]);
    double q=mxGetScalar(prhs[3]),v=0,dt=mxGetScalar(prhs[4]),count=mxGetScalar(prhs[5]);
    if(mxGetNumberOfElements(prhs[3])==2) v=mxGetPr(prhs[3])[1];
    if(!isfinite(count)||count<1||count>10000||count!=floor(count)||!isfinite(dt)||dt<=0)
        mexErrMsgIdAndTxt("hill:step","Invalid integration step.");
    int steps=(int)count; double h=dt/steps;
    plhs[0]=mxCreateDoubleMatrix(n,2,mxREAL); double *s=mxGetPr(plhs[0]);
    plhs[1]=mxCreateLogicalScalar(0); s[0]=q; s[n]=v;
    for(mwSize k=0;k<n-1;k++) {
        int failed=0;
        for(int j=0;j<steps;j++) {
            double f=(double)j/steps,da=a[k+1]-a[k],dm=m[k+1]-m[k];
            double acc1=acceleration(p,q,v,a[k]+da*f,m[k]+dm*f);
            double v2=v+h*acc1/2,acc2=acceleration(p,q+h*v/2,v2,a[k]+da*(f+.5/steps),m[k]+dm*(f+.5/steps));
            double v3=v+h*acc2/2,acc3=acceleration(p,q+h*v2/2,v3,a[k]+da*(f+.5/steps),m[k]+dm*(f+.5/steps));
            double v4=v+h*acc3,acc4=acceleration(p,q+h*v3,v4,a[k]+da*(f+1.0/steps),m[k]+dm*(f+1.0/steps));
            q+=h*(v+2*v2+2*v3+v4)/6; v+=h*(acc1+2*acc2+2*acc3+acc4)/6;
            if(!isfinite(q+v)||fabs(q)>4||fabs(v)>50) {failed=1;break;}
        }
        if(failed) {
            *((mxLogical*)mxGetData(plhs[1]))=1;
            for(mwSize kk=k+1;kk<n;kk++) s[kk]=4;
            break;
        }
        s[k+1]=q; s[n+k+1]=v;
    }
}
