/* Development-only tuning; MATLAB is the delivered implementation. */
#include "../../Demo/BSE/private/native/hill_plant.h"
static double motor(double previous,double requested,double increment,double cap) {
  requested=fmax(-cap,fmin(cap,requested));
  if(previous*requested<0) previous=0;
  return copysign(fmin(fabs(requested),fabs(previous)+increment),requested);
}
API int dual_pid(const double *p,const double *r,const double *cfg,const double *a,
    const double *mass,const double *wait,const double *threshold,const double *initial,
    const double *qref,const double *vref,int count,double dt,int steps,int stride,
    double cap,double rate,const double *g,double *out) {
  double s[STATE_SIZE],m[5]={0},command=0,period=stride*dt;
  hill_mpc_init(p,r,cfg,a,mass[0],wait[0],initial,s,out);
  for(int k=1;k<count;k++) {
    if((k-1)%stride==0) {
      double *actual=out+(k-1)*OUTPUT_SIZE,error=actual[3]-threshold[k-1];
      double d=exp(-period/g[3]); if(m[3]==0)m[3]=actual[3];
      m[2]=d*m[2]+(1-d)*(actual[3]-m[3])/period; m[3]=actual[3];
      if(wait[k-1]>.5){m[0]=0;m[2]=0;m[4]=0;}
      if(error>0 && wait[k-1]<=.5)m[4]=1;
      double qe=qref[k-1]-actual[0],ve=vref[k-1]-actual[1];
      double i0=m[0],i1=fmax(-30,fmin(30,m[1]+g[5]*period*qe)),torque=0;
      if(m[4]>0) {
        i0=fmax(0,fmin(cap,m[0]+g[1]*period*error));
        torque=fmax(0,g[0]*error+i0+g[2]*m[2]);
      }
      double raw=torque+g[4]*qe+i1+g[6]*ve;
      double reachable=motor(s[28],raw,rate*period,cap);
      if((raw-reachable)*(i0+i1-m[0]-m[1])<=0){m[0]=i0;m[1]=i1;}
      if(m[4]>0)torque=fmax(0,g[0]*error+m[0]+g[2]*m[2]);
      command=fmax(-cap,fmin(cap,torque+g[4]*qe+m[1]+g[6]*ve));
    }
    double u=motor(s[28],command,rate*dt,cap);
    if(s[28]*u<0)s[28]=0; else if(fabs(u)<fabs(s[28]))s[28]=u;
    s[7]=0;
    if(hill_mpc_advance(p,r,cfg,s,a+2*k,mass[k],wait[k],u,dt,steps,out+k*OUTPUT_SIZE))return k+1;
    s[7]=0;out[k*OUTPUT_SIZE+22]=0;
  }
  return 0;
}
