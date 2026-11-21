/* Native replay of the same reactive PID for subject/source-wide calibration. */
#include "../../Demo/NMPC/private/native/current_mtu_control.c"
#if defined(_WIN32)
#define API __declspec(dllexport)
#else
#define API
#endif
API int current_pid_run(const double *p,const double *cfg,const double *pd,const double *actuator,
                       const double *pid,int n,const double *inputs,const double *mass,const double *waiting,
                       const double *initial,const double *threshold,double dt,int steps,double *result) {
    double s[CS],row[NO];controlled_init(p,cfg,inputs,mass[0],waiting[0],initial,s,row);
    memcpy(result,row,sizeof(row));
    double integral=0,derivative=0,previous=row[3],command=0,period=2*dt;int armed=0;
    for(int k=1;k<n;k++) {
        if((k-1)%2==0) {
            double previousCommand=command;
            double measured=row[3],error=measured-threshold[k-1],decay=exp(-period/pid[3]);
            derivative=decay*derivative+(1-decay)*(measured-previous)/period;previous=measured;
            if(waiting[k-1]>.5) {command=0;integral=0;derivative=0;armed=0;}
            else {
                if(error>0)armed=1;
                if(armed) {
                    double candidate=integral+pid[1]*period*error;
                    double raw=pid[0]*error+candidate+pid[2]*derivative;
                    double reachable=fmin(actuator[0],s[28]+actuator[1]*period);
                    if(!((raw>reachable && error>0)||(raw<0 && error<0)))integral=candidate;
                    integral=fmax(0,fmin(actuator[0],integral));
                    command=fmax(0,fmin(actuator[0],pid[0]*error+integral+pid[2]*derivative));
                    if(command<1e-8 && error<=0 && integral<1e-6)armed=0;
                } else command=0;
            }
            if(waiting[k-1]<=.5 && pid[4]>0 && command<previousCommand) {
                command=command+(previousCommand-command)*exp(-period/pid[4]);
                if(command<.05)command=0;
            }
        }
        double applied=fmax(0,fmin(fmin(actuator[0],s[28]+actuator[1]*dt),command));
        if(controlled_step(p,cfg,pd,actuator,s,inputs+2*k,mass[k],waiting[k],applied,dt,steps,row))return k;
        memcpy(result+NO*k,row,sizeof(row));
    }
    return n;
}
