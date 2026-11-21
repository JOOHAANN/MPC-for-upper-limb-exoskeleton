/* Static force inversion is used ONLY to request neural commands.
 * Actual forces always come from the compliant Hill solve below. */
static double intent_input(double back, double front, double mass,
                           double difference, const double *r) {
  double load = mass / 15.0;
  return back + r[10] * front + r[11] * front * load + r[12] * back * load +
         r[13] * back * back + r[14] * difference * load;
}

static double excitation_for_torque(const double *unit, const double *p,
                                    int index, double q, double torque) {
  double lo = p[20 + index] / unit[1], w = sin(p[27]);
  double slack = p[22 + index] - lo * sqrt(unit[2] * unit[2] - w * w);
  double lmt =
      p[22 + index] - (index == 0 ? 1 : -1) * p[20 + index] * (q - 1.4);
  double ft = fmax(0, torque / unit[0]);
  double lt = slack * (1 + unit[6] * sqrt(ft));
  double x = (lmt - lt) / lo;
  if (x <= 0)
    return 0; /* No feasible tensile force at this slack geometry. */
  double l = sqrt(x * x + w * w), cs = x / l, fl = exp(-pow((l - 1) / .45, 2));
  double activation = (ft / cs - passive(l)) / fmax(fl, 1e-8);
  return pow(clamp01(activation), 1 / unit[5]);
}

static void stance_command(const double *p, const double *front,
                           const double *cfg, double q, double v, double mass,
                           double integral, double *command, double *desired,
                           double *correction) {
  double passive_torque[2];
  for (int u = 0; u < 2; u++) {
    double x, ft, ce, pe, vn, res;
    const double *unit = u == 0 ? p : front;
    fiber(unit, q, 0, 0, .001, 1, 0, u == 0 ? 1 : -1, p[20 + u], p[22 + u],
          sin(p[27]), &x, &ft, &ce, &pe, &vn, &res);
    passive_torque[u] = unit[0] * ft;
  }
  *correction =
      cfg[4] *
      tanh((cfg[1] * (cfg[0] - q) - cfg[2] * v + cfg[3] * integral) / cfg[4]);
  *desired = -gravity(p, q, mass) + *correction;
  double difference = *desired - passive_torque[0] + passive_torque[1];
  double target[2] = {passive_torque[0] + cfg[5] + fmax(difference, 0),
                      passive_torque[1] + cfg[5] + fmax(-difference, 0)};
  command[0] = excitation_for_torque(p, p, 0, q, target[0]);
  command[1] = excitation_for_torque(front, p, 1, q, target[1]);
}

static double decode_intent_signal(double a, double lead, double tonic,
                                   double gain) {
  return clamp01(a + lead - gain * tonic);
}
static void movement_command(const double *p, const double *front,
                             const double *r, double q, double v, double mass,
                             double intent, double front_input, double *command,
                             double *desired, double *correction) {
  double minimum[2], base_front;
  for (int u = 0; u < 2; u++) {
    double x, ft, ce, pe, vn, res;
    const double *unit = u == 0 ? p : front;
    fiber(unit, q, 0, 0, .001, 1, 0, u == 0 ? 1 : -1, p[20 + u], p[22 + u],
          sin(p[27]), &x, &ft, &ce, &pe, &vn, &res);
    minimum[u] = unit[0] * ft;
  }
  double x, ft, ce, pe, vn, res;
  fiber(front, q, pow(clamp01(front_input), front[5]), 0, .001, 1, 0, -1, p[21],
        p[23], sin(p[27]), &x, &ft, &ce, &pe, &vn, &res);
  base_front = front[0] * ft;
  *correction = r[3] * (intent - q) - r[4] * v;
  double request = -gravity(p, q, mass) + *correction;
  /* Asymmetric request saturation permits gravity-driven flexion with a
   * modest antagonist request. It does not cap actual muscle torque. */
  *desired = request < 0 ? r[5] * tanh(request / r[5]) : request;
  double back_target = fmax(minimum[0], *desired + base_front);
  double front_target = fmax(base_front, minimum[0] - *desired);
  back_target += r[7];
  front_target += r[7];
  command[0] = excitation_for_torque(p, p, 0, q, back_target);
  command[1] = excitation_for_torque(front, p, 1, q, front_target);
}
