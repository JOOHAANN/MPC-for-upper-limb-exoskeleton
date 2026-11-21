/* Development-only PID tuning. MATLAB remains the delivered runtime. */
#include "../../Demo/BSE/private/native/hill_plant.h"

/* gains: Kp, Ki, Kd, derivative time constant, integral release time. */
API int compare_pid_rollout(const double *p, const double *r, const double *cfg,
                           const double *activation, const double *mass,
                           const double *waiting, const double *threshold,
                           const double *initial, int count, double dt,
                           int substeps, int stride, double cap, double rate,
                           const double *gains, double *out) {
    double state[STATE_SIZE], integral = 0, derivative = 0, command = 0;
    int active = 0;
    if (hill_mpc_init(p,r,cfg,activation,mass[0],waiting[0],initial,state,out)) return 1;
    double previous = out[3], period = stride*dt;
    for (int k=1;k<count;k++) {
        if ((k-1)%stride == 0) {
            double measured = out[(k-1)*OUTPUT_SIZE+3];
            double error = measured-threshold[k-1];
            double decay = exp(-period/gains[3]);
            derivative = decay*derivative + (1-decay)*(measured-previous)/period;
            previous = measured;
            if (waiting[k-1] > .5) {
                command=0; integral=0; derivative=0; active=0;
            } else {
                if (error>0) active=1;
                if (active) {
                    double trial = integral+gains[1]*period*error;
                    double raw = gains[0]*error+trial+gains[2]*derivative;
                    if (!((raw>cap && error>0) || (raw<0 && error<0))) integral=trial;
                    integral=fmax(0,fmin(cap,integral));
                    command=fmax(0,fmin(cap,gains[0]*error+integral+gains[2]*derivative));
                    if (command<1e-8 && error<=0 && integral<1e-6) active=0;
                } else command=0;
            }
        }
        double applied=fmax(fmax(0,state[28]-rate*dt),
                            fmin(fmin(cap,state[28]+rate*dt),command));
        if (hill_mpc_advance(p,r,cfg,state,activation+2*k,mass[k],waiting[k],
                             applied,dt,substeps,out+k*OUTPUT_SIZE)) return k+1;
    }
    return 0;
}
