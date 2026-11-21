/* Compliant Hill muscle coupled to joint forward dynamics. No angle target.
 * Backward Euler is solved in the fiber and joint at every integration step.
 * Tendon force is monotone and tensile; no polynomial stiffness or force cap. */
#include <math.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
static double clamp01(double x)
{
    return fmax(0, fmin(1, x));
}
static double fv(double v)
{
    return v <= 0 ? fmax(0, 1 + v) / (1 - v / .25) : (1 + 1.5 * v / .25) / (1 + v / .25);
}
static double passive(double l)
{
    return fmax(0, expm1(fmin(35, 5 * (l - 1) / .6)) / expm1(5));
}
static void fiber(const double *p, double q, double a, double oldl, double h, int initial, double v0, int direction,
                  double R, double Lref, double w, double *xn, double *ft, double *ce, double *pe, double *vn,
                  double *res)
{
    (void)v0; /* Joint speed is independent of initial compliant-fiber speed. */
    double lo = R / p[1], xr = sqrt(p[2] * p[2] - w * w);
    double slack = Lref - lo * xr, c = p[6] * slack / lo, ell = xr - direction * p[1] * (q - 1.4);
    /* Zero projected fiber length is the slack-tendon geometric limit.
     * A positive numerical floor would leave a spurious transmitted force. */
    /* A slack tendon may shorten below its slack length. The upper fiber
     * bound is the full available MT length, not only one tendon strain.
     * The previous ell+c bracket falsely excluded valid slack-tendon roots. */
    double lower = 0, upper = fmax(1e-5, ell + slack / lo),
           x = initial ? (lower + upper) / 2 : sqrt(fmax(0, oldl * oldl - w * w));
    x = fmax(lower + 1e-10, fmin(upper - 1e-10, x));
    for (int j = 0; j < 40; j++)
    {
        double l = sqrt(x * x + w * w), cs = x / l, dcs = w * w / (l * l * l);
        /* A compliant tendon can change length while the initial fiber is
         * stationary. Initialize fiber equilibrium at zero fiber velocity;
         * the original joint velocity remains an independent initial state. */
        double vel = initial ? 0 : (l - oldl) / (h * p[4]);
        double dvel = initial ? 0 : cs / (h * p[4]);
        double stretch = fmax(0, ell - x), tendon = stretch * stretch / (c * c), dtendon = -2 * stretch / (c * c);
        double fl = exp(-pow((l - 1) / .45, 2)), dfl = -2 * (l - 1) / (.45 * .45) * fl * cs;
        double fvel = fv(vel), dfvel = vel <= -1 ? 0 : (vel <= 0 ? 5 / pow(1 - 4 * vel, 2) : 2 / pow(1 + 4 * vel, 2));
        double fp = passive(l),
               dfp = l <= 1 || 5 * (l - 1) / .6 >= 35 ? 0 : (5 / .6) * exp(5 * (l - 1) / .6) * cs / expm1(5);
        double active = a * fl * fvel, dactive = a * (dfl * fvel + fl * dfvel * dvel);
        double rhs = (active + fp + .03 * vel) * cs;
        double f = tendon - rhs, der = dtendon - (dactive + dfp + .03 * dvel) * cs - (active + fp + .03 * vel) * dcs;
        if (fabs(f) < 3e-11)
            break;
        if (f > 0)
            lower = x;
        else
            upper = x;
        double next = x - f / der;
        if (!isfinite(next) || next <= lower || next >= upper)
            next = (lower + upper) / 2;
        x = next;
    }
    double l = sqrt(x * x + w * w), cs = x / l;
    *vn = initial ? 0 : (l - oldl) / (h * p[4]);
    *ft = pow(fmax(0, ell - x) / c, 2);
    *ce = a * exp(-pow((l - 1) / .45, 2)) * fv(*vn) * cs;
    *pe = passive(l) * cs;
    *xn = x;
    *res = *ft - *ce - *pe - .03 * (*vn) * cs;
}
static double gravity(const double *p, double q, double m)
{
    return -9.81 * p[18] * p[19] * cos(q) - 9.81 * m * p[12] * cos(q - p[13]);
}
/* Elastic dead-zone soft limit. C1 onset, no angle clipping, no reference
 * feedback and no boundary damping. Its torque is NOT part of Hill torque. */
static double boundary(const double *p, double q)
{
    double low = fmax(p[30] - q, 0), high = fmax(q - p[31], 0);
    low = low * low / (low + .05);
    high = high * high / (high + .05);
    return p[29] * (tanh(p[28] * low / p[29]) - tanh(p[28] * high / p[29]));
}
