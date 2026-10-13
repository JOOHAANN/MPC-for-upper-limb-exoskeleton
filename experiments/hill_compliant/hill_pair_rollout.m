function p = hill_pair_rollout(par, inputs, m, initial, dt, steps, useMex)
    % HILL_PAIR_ROLLOUT Opposing compliant Hill muscles drive one joint.
    % Inputs are back/front activations. No reference angle/torque is accepted.
    if nargin < 6
        steps = 2;
    end
    if nargin < 7
        useMex = true;
    end
    par = par(:)';
    assert(numel(par) == 32 && all(isfinite(par)) && all(par([1 2 3 5 6 7 8 9 10 11 12 13 19:30]) > 0));
    assert(par(31) < par(32) && par(28) > 0 && par(28) < pi / 2);
    assert(all(par(23:24) - par(21:22) ./ par([2 11]) .* sqrt(par([3 12]).^2 - sin(par(28))^2) > 0));
    assert(size(inputs, 2) == 2 && size(inputs, 1) == numel(m) && all(isfinite(inputs(:))));
    assert(iscolumn(m) && all(m >= 0) && numel(initial) == 2 && all(isfinite(initial)));
    assert(dt > 0 && steps >= 1 && steps == floor(steps));
    if useMex && exist('hill_pair_mex', 'file') == 3
        [x, failed] = hill_pair_mex(par, inputs, m, initial, dt, steps);
    else
        [x, failed] = matlab_integrate(par, inputs, m, initial, dt, steps);
    end
    q = x(:, 1);
    v = x(:, 2);
    gravity = -9.81 * par(19) * par(20) * cos(q) - 9.81 * m * par(13) .* cos(q - par(14));
    I = m * par(13)^2 + 4 * par(19) * par(20)^2 / 3;
    boundary = hill_soft_limit(q, par);
    p = struct('q', q, 'v', v, 'torque', x(:, 3), 'backTorque', x(:, 4), 'frontTorque', x(:, 5), ...
               'backActive', x(:, 6), 'backPassive', x(:, 7), 'backFiberDamping', x(:, 8), ...
               'frontActive', x(:, 9), 'frontPassive', x(:, 10), 'frontFiberDamping', x(:, 11), ...
               'backLengthNormalized', x(:, 12), 'frontLengthNormalized', x(:, 13), ...
               'effectiveActivation', x(:, 14:15), 'equilibriumResidual', x(:, 16), ...
               'momentumResidual', x(:, 17), 'gravityTorque', gravity, 'dampingTorque', -par(4) * v, ...
               'acceleration', (x(:, 3) + gravity - par(4) * v + boundary) ./ I, 'boundaryTorque', boundary, ...
               'failed', failed);
end

function [out, failed] = matlab_integrate(par, a, m, initial, dt, steps)
    front = [par(9), par(11), par(12), 0, par(25), par(10), par(26), par(27)];
    units = [par(1:8); front];
    q = initial(1);
    v = initial(2);
    h = dt / steps;
    n = size(a, 1);
    w = sin(par(28));
    gamma = reshape(units(:, 6), 1, 2);
    activationTau = reshape(units(:, 8), 1, 2);
    act = max(0, min(1, a(1, :))).^gamma;
    fiberState = zeros(2, 6);
    oldl = zeros(1, 2);
    for u = 1:2
        fiberState(u, :) = fiber(units(u, :), q, act(u), 0, h, true, v, 3 - 2 * u, par(20 + u), par(22 + u), w);
        oldl(u) = sqrt(fiberState(u, 1)^2 + w^2);
    end
    out = zeros(n, 17);
    failed = false;
    for k = 1:n
        mom = 0;
        if k > 1
            for s = 1:steps
                f = s / steps;
                mass = m(k - 1) + f * (m(k) - m(k - 1));
                target = max(0, min(1, a(k - 1, :) + f * (a(k, :) - a(k - 1, :)))).^gamma;
                act = target + (act - target) .* exp(-h ./ activationTau);
                I = mass * par(13)^2 + 4 * par(19) * par(20)^2 / 3;
                qold = q;
                vold = v;
                guess = q + h * v;
                for j = 1:9
                    for u = 1:2
                        fiberState(u, :) = fiber(units(u, :), guess, act(u), oldl(u), h, false, vold, 3 - 2 * u, par(20 + u), par(22 + u), w);
                    end
                    muscle = par(1) * fiberState(1, 2) - par(9) * fiberState(2, 2);
                    fun = I * (guess - qold - h * vold) / h^2 + par(4) * (guess - qold) / h - muscle - gravity(guess, mass, par) - hill_soft_limit(guess, par);
                    if abs(fun) < 1e-7
                        break
                    end
                    ep = 1e-6;
                    tplus = 0;
                    for u = 1:2
                        fs = fiber(units(u, :), guess + ep, act(u), oldl(u), h, false, vold, 3 - 2 * u, par(20 + u), par(22 + u), w);
                        tplus = tplus + (3 - 2 * u) * units(u, 1) * fs(2);
                    end
                    derivative = I / h^2 + par(4) / h - (tplus - muscle) / ep - (gravity(guess + ep, mass, par) - gravity(guess, mass, par)) / ep ...
                                 - (hill_soft_limit(guess + ep, par) - hill_soft_limit(guess, par)) / ep;
                    guess = guess - max(-.1, min(.1, fun / derivative));
                end
                q = guess;
                v = (q - qold) / h;
                for u = 1:2
                    fiberState(u, :) = fiber(units(u, :), q, act(u), oldl(u), h, false, vold, 3 - 2 * u, par(20 + u), par(22 + u), w);
                    oldl(u) = sqrt(fiberState(u, 1)^2 + w^2);
                end
                mom = I * (v - vold) / h - par(1) * fiberState(1, 2) + par(9) * fiberState(2, 2) - gravity(q, mass, par) - hill_soft_limit(q, par) + par(4) * v;
                if ~isfinite(q + v) || abs(q) > 4 || abs(v) > 50 || sum(abs(fiberState(:, 6))) > 1e-4 || abs(mom) > 1
                    failed = true;
                    break
                end
            end
        end
        force = fiberState(:, 2) .* units(:, 1);
        ce = fiberState(:, 3) .* units(:, 1);
        pe = fiberState(:, 4) .* units(:, 1);
        damp = units(:, 1) * .03 .* fiberState(:, 5) .* fiberState(:, 1) ./ oldl';
        out(k, :) = [q, v, force(1) - force(2), force', ce(1), pe(1), damp(1), ce(2), pe(2), damp(2), ...
                     oldl, act, max(abs(fiberState(:, 6))), mom];
        if failed
            out(k + 1:end, 1) = 4;
            break
        end
    end
end

function fs = fiber(p, q, a, oldl, h, initial, v0, direction, R, Lref, w)
    lo = R / p(2);
    xr = sqrt(p(3)^2 - w * w);
    slack = Lref - lo * xr;
    c = p(7) * slack / lo;
    ell = xr - direction * p(2) * (q - 1.4);
    lower = 0;
    % A slack tendon can be shorter than its slack length. Allow the fiber
    % projection to span the full available MT length, with nonnegative LT.
    upper = max(1e-5, ell + slack / lo);
    for j = 1:36
        x = (lower + upper) / 2;
        l = sqrt(x * x + w * w);
        cs = x / l;
        if initial
            % Joint velocity is preserved; a compliant/slack tendon allows
            % independent fiber velocity at initialization.
            vn = 0;
        else
            vn = (l - oldl) / (h * p(5));
        end
        ft = (max(0, ell - x) / c)^2;
        ce = a * exp(-((l - 1) / .45)^2) * force_velocity(vn);
        rhs = (ce + passive(l) + .03 * vn) * cs;
        if ft > rhs
            lower = x;
        else
            upper = x;
        end
    end
    l = sqrt(x * x + w * w);
    cs = x / l;
    if initial
        vn = 0;
    else
        vn = (l - oldl) / (h * p(5));
    end
    ft = (max(0, ell - x) / c)^2;
    ce = a * exp(-((l - 1) / .45)^2) * force_velocity(vn) * cs;
    pe = passive(l) * cs;
    fs = [x, ft, ce, pe, vn, ft - ce - pe - .03 * vn * cs];
end

function y = force_velocity(v)
    if v <= 0
        y = max(0, 1 + v) / (1 - v / .25);
    else
        y = (1 + 1.5 * v / .25) / (1 + v / .25);
    end
end

function y = passive(l)
    y = max(0, expm1(min(35, 5 * (l - 1) / .6)) / expm1(5));
end

function g = gravity(q, m, p)
    g = -9.81 * p(19) * p(20) * cos(q) - 9.81 * m * p(13) * cos(q - p(14));
end
