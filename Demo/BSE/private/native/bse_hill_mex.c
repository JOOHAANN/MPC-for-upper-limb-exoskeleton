/* MATLAB gateway to the shared Hill plant. No references enter this file. */
#include "hill_plant.h"
#include "mex.h"

static void require_vector(const mxArray *a, mwSize n, const char *name) {
  if (!mxIsDouble(a) || mxIsComplex(a) || mxIsSparse(a) ||
      mxGetNumberOfElements(a) != n)
    mexErrMsgIdAndTxt("BSE:Input", "%s must contain %u real doubles.", name,
                      (unsigned)n);
  const double *x = mxGetPr(a);
  for (mwSize k = 0; k < n; ++k)
    if (!mxIsFinite(x[k]))
      mexErrMsgIdAndTxt("BSE:Input", "%s must be finite.", name);
}

static void copy_row(double *target, const double *row, mwSize n, mwSize k) {
  for (mwSize j = 0; j < OUTPUT_SIZE; ++j)
    target[k + n * j] = row[j];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  char action[24];
  if (nrhs < 4 || !mxIsChar(prhs[0]) ||
      mxGetString(prhs[0], action, sizeof(action)))
    mexErrMsgIdAndTxt("BSE:Action", "Choose init, step, preview or rollout.");
  require_vector(prhs[1], 32, "Physical parameters");
  mwSize neuralCount = mxGetNumberOfElements(prhs[2]);
  if (neuralCount != 10 && neuralCount != 15 && neuralCount != 16)
    mexErrMsgIdAndTxt("BSE:Input", "Expected 10, 15 or 16 neural parameters.");
  require_vector(prhs[2], neuralCount, "Neural parameters");
  require_vector(prhs[3], 7, "Standing parameters");
  double neural[16] = {0};
  memcpy(neural, mxGetPr(prhs[2]), neuralCount * sizeof(double));
  const double *p = mxGetPr(prhs[1]), *r = neural, *cfg = mxGetPr(prhs[3]);
  int failed = 0;

  if (!strcmp(action, "init")) {
    if (nrhs != 8 || nlhs != 2)
      mexErrMsgIdAndTxt("BSE:Arguments", "Invalid init call.");
    require_vector(prhs[4], 2, "Activation");
    require_vector(prhs[5], 1, "Mass");
    require_vector(prhs[6], 1, "Waiting");
    require_vector(prhs[7], 2, "Initial state");
    plhs[0] = mxCreateDoubleMatrix(STATE_SIZE, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, OUTPUT_SIZE, mxREAL);
    failed = hill_mpc_init(p, r, cfg, mxGetPr(prhs[4]), mxGetScalar(prhs[5]),
                           mxGetScalar(prhs[6]), mxGetPr(prhs[7]),
                           mxGetPr(plhs[0]), mxGetPr(plhs[1]));
  } else if (!strcmp(action, "step")) {
    if (nrhs != 11 || nlhs != 2)
      mexErrMsgIdAndTxt("BSE:Arguments", "Invalid step call.");
    require_vector(prhs[4], STATE_SIZE, "Plant state");
    require_vector(prhs[5], 2, "Activation");
    for (int j = 6; j < 11; ++j)
      require_vector(prhs[j], 1, "Step argument");
    double dt = mxGetScalar(prhs[9]);
    int steps = (int)mxGetScalar(prhs[10]);
    if (dt <= 0 || steps < 1)
      mexErrMsgIdAndTxt("BSE:Step", "Positive step sizes required.");
    plhs[0] = mxDuplicateArray(prhs[4]);
    plhs[1] = mxCreateDoubleMatrix(1, OUTPUT_SIZE, mxREAL);
    failed =
        hill_mpc_advance(p, r, cfg, mxGetPr(plhs[0]), mxGetPr(prhs[5]),
                         mxGetScalar(prhs[6]), mxGetScalar(prhs[7]),
                         mxGetScalar(prhs[8]), dt, steps, mxGetPr(plhs[1]));
  } else if (!strcmp(action, "rollout")) {
    if (nrhs != 10 || nlhs != 1)
      mexErrMsgIdAndTxt("BSE:Arguments", "Invalid rollout call.");
    mwSize n = mxGetM(prhs[4]);
    if (n < 1 || mxGetN(prhs[4]) != 2)
      mexErrMsgIdAndTxt("BSE:Activation", "Activation must be N-by-2.");
    require_vector(prhs[4], 2 * n, "Activation");
    require_vector(prhs[5], n, "Mass");
    require_vector(prhs[6], n, "Waiting");
    require_vector(prhs[7], 2, "Initial state");
    require_vector(prhs[8], 1, "Sample period");
    require_vector(prhs[9], 1, "Substeps");
    double dt = mxGetScalar(prhs[8]);
    int steps = (int)mxGetScalar(prhs[9]);
    if (dt <= 0 || steps < 1)
      mexErrMsgIdAndTxt("BSE:Step", "Positive step sizes required.");
    const double *a = mxGetPr(prhs[4]), *mass = mxGetPr(prhs[5]),
                 *wait = mxGetPr(prhs[6]);
    double state[STATE_SIZE], row[OUTPUT_SIZE], input[2] = {a[0], a[n]};
    plhs[0] = mxCreateDoubleMatrix(n, OUTPUT_SIZE, mxREAL);
    hill_mpc_init(p, r, cfg, input, mass[0], wait[0], mxGetPr(prhs[7]), state,
                  row);
    copy_row(mxGetPr(plhs[0]), row, n, 0);
    for (mwSize k = 1; k < n; ++k) {
      input[0] = a[k];
      input[1] = a[n + k];
      failed = hill_mpc_advance(p, r, cfg, state, input, mass[k], wait[k], 0,
                                dt, steps, row);
      if (failed)
        break;
      copy_row(mxGetPr(plhs[0]), row, n, k);
    }
  } else if (!strcmp(action, "preview")) {
    if (nrhs != 12 || nlhs != 1)
      mexErrMsgIdAndTxt("BSE:Arguments", "Invalid preview call.");
    require_vector(prhs[4], STATE_SIZE, "Plant state");
    mwSize blocks = mxGetNumberOfElements(prhs[5]), n = mxGetM(prhs[6]);
    if (!blocks || !n || mxGetN(prhs[6]) != 2)
      mexErrMsgIdAndTxt("BSE:Preview", "Invalid prediction dimensions.");
    require_vector(prhs[5], blocks, "Control blocks");
    require_vector(prhs[6], 2 * n, "Activation forecast");
    require_vector(prhs[7], n, "Mass forecast");
    require_vector(prhs[8], n, "Waiting forecast");
    for (int j = 9; j < 12; ++j)
      require_vector(prhs[j], 1, "Prediction argument");
    double dt = mxGetScalar(prhs[9]), rate = mxGetScalar(prhs[11]);
    int steps = (int)mxGetScalar(prhs[10]);
    if (dt <= 0 || steps < 1 || rate <= 0)
      mexErrMsgIdAndTxt("BSE:Step", "Positive prediction settings required.");
    double *activation = mxCalloc(2 * n, sizeof(double));
    double *out = mxCalloc(n * OUTPUT_SIZE, sizeof(double));
    const double *a = mxGetPr(prhs[6]);
    for (mwSize k = 0; k < n; ++k) {
      activation[2 * k] = a[k];
      activation[2 * k + 1] = a[n + k];
    }
    failed =
        hill_mpc_predict_map(p, r, cfg, mxGetPr(prhs[4]), mxGetPr(prhs[5]),
                             activation, mxGetPr(prhs[7]), mxGetPr(prhs[8]),
                             (int)blocks, (int)n, dt, steps, rate, out);
    plhs[0] = mxCreateDoubleMatrix(failed ? 0 : n, OUTPUT_SIZE, mxREAL);
    if (!failed)
      for (mwSize k = 0; k < n; ++k)
        copy_row(mxGetPr(plhs[0]), out + k * OUTPUT_SIZE, n, k);
    mxFree(activation);
    mxFree(out);
    return;
  } else {
    mexErrMsgIdAndTxt("BSE:Action", "Unknown plant operation.");
  }
  if (failed)
    mexErrMsgIdAndTxt("BSE:Equilibrium",
                      "The coupled Hill solve did not converge.");
}
