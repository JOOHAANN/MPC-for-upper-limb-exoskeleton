"""Create the controlled 24-parameter core from the verified BSE implementation.

The shared Hill mechanics header is included directly from the current BSE.
Only assistance and explicitly gated activation feedback are added here.
"""
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
source=ROOT/'Demo/BSE/Direct_Activation_No_Limits/private/native/direct_mtu_mex.c'
text=source.read_text();text=text[:text.index('#if defined(_WIN32)')]
text=text.replace('#include "mtu_mechanics.h"','#include "../../../BSE/Direct_Activation_No_Limits/private/native/mtu_mechanics.h"')
text=text.replace('#define NS 34','#define NS 34\n#define CS 69')
text=text.replace('double *command,double *request,double *correction)',
                  'double assist,double *command,double *request,double *correction)')
text=text.replace('*request=-gravity(p,q,mass)+*correction;',
                  '*request=-gravity(p,q,mass)+*correction-assist;')
text=text.replace('mass,0,cmd,&s[11],&s[12])','mass,0,0,cmd,&s[11],&s[12])')
text=text.replace('static int advance(', 'static int advance_physical(')
text=text.replace('double massEnd,double waitEnd,double dt,int steps,double *o)',
    'double massEnd,double waitEnd,double assistEnd,const double *target,\n'+
    '                   const double *pd,double *filtered,double dt,int steps,double *o)')
text=text.replace('double h=dt/steps,w=sin(p[23]);','double h=dt/steps,w=sin(p[23]);\n'+
    '    /* Release is instantaneous; an increase follows the external rise limit. */\n'+
    '    double startAssist=assistEnd<s[28]?assistEnd:s[28];')
text=text.replace('double f=(double)k/steps,mass=s[32]+f*(massEnd-s[32]);',
    'double f=(double)k/steps,mass=s[32]+f*(massEnd-s[32]);\n'+
    '        double assist=startAssist+f*(assistEnd-startAssist);')
text=text.replace('mass,s[7],cmd,&request,&correction)', 'mass,s[7],assist,cmd,&request,&correction)')
needle='        s[11]=s[6]*request;s[12]=s[6]*correction;'
replacement='''        /* Original movement PD: correct the BACK activation command only.
         * Targets are current causal unassisted-model states, never IK data.
         * The current measured assistance, not the requested command, gates PD. */
        if(waitEnd<=.5 && fabs(assist)>1e-9 && target && pd) {
            double delta=pd[0]*(target[2*(k-1)]-s[0])+pd[1]*(target[2*(k-1)+1]-s[1]);
            if(s[0]<1)delta=fmax(-pd[3],fmin(pd[2],delta));
            *filtered=pd[4]>0?delta+(*filtered-delta)*exp(-h/pd[4]):delta;
            input[0]=clamp01(input[0]+*filtered);
        } else *filtered=0;
'''+needle
assert needle in text;text=text.replace(needle,replacement)
text=text.replace('-muscle-gravity(p,guess,mass);','-muscle-gravity(p,guess,mass)-assist;')
text=text.replace('-gravity(p,s[0],mass)+p[3]*s[1];','-gravity(p,s[0],mass)+p[3]*s[1]-assist;')
text=text.replace('s[30]=a[0];s[31]=a[1];s[32]=massEnd;s[33]=waitEnd;',
                  's[28]=assistEnd;s[30]=a[0];s[31]=a[1];s[32]=massEnd;s[33]=waitEnd;')
text=text.replace('/* Columns 24/25/26 remain zero: no intent, assistance or angle boundary. */',
                  '/* The controlled wrapper fills target/assistance columns. Boundary is zero. */')
text=text.replace('/* Standing can blend in smoothly. It exits immediately at an action\n         * sample: zero posture feedback is guaranteed throughout movement. */',
                  '/* Standing feedback blends in and exits at the first action sample.\n         * Movement activation PD is handled separately below. */')
text=text.replace('/* Direct activation -> two compliant MTUs -> one joint.\n * Standing balance is the only human feedback. No movement intent, movement\n * PD, angle boundary torque, angle clipping or reference trajectory exists. */',
'''/* Current calibrated dual-MTU model with external assistance.
 * At zero assistance, this is the verified BSE with standing feedback only.
 * Movement PD is gated by applied torque and corrects back activation.
 * Physical coefficients, Hill mechanics and angle-limit removal are unchanged. */''')
text+='''
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
'''
target=ROOT/'Demo/NMPC/private/native/current_mtu_control.c'
target.write_text(text)
print('Wrote',target)
