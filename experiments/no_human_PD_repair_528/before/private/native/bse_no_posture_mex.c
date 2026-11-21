/* Separate ablation gateway. The shared Hill equations are never edited. */
#include "../../../BSE/private/native/hill_plant.h"
#include "mex.h"

static void vector(const mxArray *a,mwSize n) {
  if (!mxIsDouble(a) || mxIsComplex(a) || mxIsSparse(a) || mxGetNumberOfElements(a)!=n)
    mexErrMsgIdAndTxt("NoPosture:Input","Invalid real-double input dimensions.");
  for(mwSize j=0;j<n;j++) if(!mxIsFinite(mxGetPr(a)[j]))
    mexErrMsgIdAndTxt("NoPosture:Input","Inputs must be finite.");
}
static void row_copy(double *out,const double *row,mwSize n,mwSize k) {
  for(mwSize j=0;j<OUTPUT_SIZE;j++) out[k+n*j]=row[j];
}
static int advance(const double *p,const double *r,const double *cfg,double *s,
                   const double *a,double mass,double wait,double u,double dt,int steps,double *out) {
  /* The inactive posture integrator is removed, rather than merely logged.
   * For signed torque, increasing magnitude is ramped; withdrawal is immediate. */
  s[7]=0;
  if (s[28]*u<0) s[28]=0;
  else if (fabs(u)<fabs(s[28])) s[28]=u;
  int failed=hill_mpc_advance(p,r,cfg,s,a,mass,wait,u,dt,steps,out);
  s[7]=0; out[22]=0;
  return failed;
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  char action[24];
  if(nrhs<4 || mxGetString(prhs[0],action,sizeof(action)))
    mexErrMsgIdAndTxt("NoPosture:Action","Choose init, step, rollout or preview.");
  vector(prhs[1],32); vector(prhs[2],16); vector(prhs[3],7);
  const double *p=mxGetPr(prhs[1]),*r=mxGetPr(prhs[2]),*cfg=mxGetPr(prhs[3]);
  if(r[3]!=0 || r[4]!=0 || cfg[1]!=0 || cfg[2]!=0 || cfg[3]!=0)
    mexErrMsgIdAndTxt("NoPosture:Feedback","Human posture gains must all be zero.");
  int failed=0;
  if(!strcmp(action,"init")) {
    if(nrhs!=8 || nlhs!=2) mexErrMsgIdAndTxt("NoPosture:Call","Invalid init call.");
    vector(prhs[4],2); vector(prhs[5],1); vector(prhs[6],1); vector(prhs[7],2);
    plhs[0]=mxCreateDoubleMatrix(STATE_SIZE,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(1,OUTPUT_SIZE,mxREAL);
    failed=hill_mpc_init(p,r,cfg,mxGetPr(prhs[4]),mxGetScalar(prhs[5]),mxGetScalar(prhs[6]),
        mxGetPr(prhs[7]),mxGetPr(plhs[0]),mxGetPr(plhs[1]));
  } else if(!strcmp(action,"step")) {
    if(nrhs!=11 || nlhs!=2) mexErrMsgIdAndTxt("NoPosture:Call","Invalid step call.");
    vector(prhs[4],STATE_SIZE); vector(prhs[5],2);
    for(int j=6;j<11;j++) vector(prhs[j],1);
    double dt=mxGetScalar(prhs[9]); int steps=(int)mxGetScalar(prhs[10]);
    if(dt<=0 || steps<1) mexErrMsgIdAndTxt("NoPosture:Step","Invalid time step.");
    plhs[0]=mxDuplicateArray(prhs[4]); plhs[1]=mxCreateDoubleMatrix(1,OUTPUT_SIZE,mxREAL);
    failed=advance(p,r,cfg,mxGetPr(plhs[0]),mxGetPr(prhs[5]),mxGetScalar(prhs[6]),
        mxGetScalar(prhs[7]),mxGetScalar(prhs[8]),dt,steps,mxGetPr(plhs[1]));
  } else if(!strcmp(action,"rollout") || !strcmp(action,"preview")) {
    int preview=!strcmp(action,"preview");
    if(nlhs!=1 || nrhs!=(preview?12:10)) mexErrMsgIdAndTxt("NoPosture:Call","Invalid rollout call.");
    int ai=preview?6:4; mwSize n=mxGetM(prhs[ai]);
    if(!n || mxGetN(prhs[ai])!=2) mexErrMsgIdAndTxt("NoPosture:Input","Excitation must be N by 2.");
    vector(prhs[ai],2*n); vector(prhs[ai+1],n); vector(prhs[ai+2],n);
    const double *a=mxGetPr(prhs[ai]),*mass=mxGetPr(prhs[ai+1]),*wait=mxGetPr(prhs[ai+2]);
    double dt=mxGetScalar(prhs[preview?9:8]); int steps=(int)mxGetScalar(prhs[preview?10:9]);
    if(dt<=0 || steps<1) mexErrMsgIdAndTxt("NoPosture:Step","Invalid time step.");
    double s[STATE_SIZE],out[OUTPUT_SIZE],input[2]={a[0],a[n]};
    mwSize blocks=0; const double *control=NULL; double rate=0;
    if(preview) {
      vector(prhs[4],STATE_SIZE); memcpy(s,mxGetPr(prhs[4]),sizeof(s));
      blocks=mxGetNumberOfElements(prhs[5]); if(!blocks) mexErrMsgIdAndTxt("NoPosture:Input","Empty controls.");
      vector(prhs[5],blocks); control=mxGetPr(prhs[5]); rate=mxGetScalar(prhs[11]);
      if(rate<=0) mexErrMsgIdAndTxt("NoPosture:Step","Positive rise rate required.");
    } else {
      vector(prhs[7],2);
      hill_mpc_init(p,r,cfg,input,mass[0],wait[0],mxGetPr(prhs[7]),s,out);
    }
    plhs[0]=mxCreateDoubleMatrix(n,OUTPUT_SIZE,mxREAL);
    if(!preview) row_copy(mxGetPr(plhs[0]),out,n,0);
    for(mwSize k=preview?0:1;k<n;k++) {
      input[0]=a[k]; input[1]=a[n+k]; double u=0;
      if(preview) {
        u=control[k*blocks/n];
        if(u*s[28]<0) s[28]=0;
        if(fabs(u)>fabs(s[28])) u=copysign(fmin(fabs(u),fabs(s[28])+rate*dt),u);
      }
      failed=advance(p,r,cfg,s,input,mass[k],wait[k],u,dt,steps,out);
      if(failed) break;
      row_copy(mxGetPr(plhs[0]),out,n,k);
    }
    if(failed && preview) {mxDestroyArray(plhs[0]);plhs[0]=mxCreateDoubleMatrix(0,OUTPUT_SIZE,mxREAL);return;}
  } else mexErrMsgIdAndTxt("NoPosture:Action","Unknown operation.");
  if(failed) mexErrMsgIdAndTxt("NoPosture:Equilibrium","The coupled Hill solve failed.");
}
