function p = hill_pair_rollout(par, inputs, m, initial, dt, steps, useMex, stance)
    % HILL_PAIR_ROLLOUT Opposing compliant Hill muscles drive one joint.
    % Inputs are back/front activations. No reference angle/torque is accepted.
    if nargin < 6
        steps = 2;
    end
    if nargin < 7
        useMex = true;
    end
    if nargin < 8
        stance = [];
    end
    par = par(:)';
    assert(numel(par) == 32 && all(isfinite(par)) && all(par([1 2 3 5 6 7 8 9 10 11 12 13 19:30]) > 0));
    assert(par(31) < par(32) && par(28) > 0 && par(28) < pi / 2);
    assert(all(par(23:24) - par(21:22) ./ par([2 11]) .* sqrt(par([3 12]).^2 - sin(par(28))^2) > 0));
    assert(size(inputs, 2) == 2 && size(inputs, 1) == numel(m) && all(isfinite(inputs(:))));
    assert(iscolumn(m) && all(m >= 0) && numel(initial) == 2 && all(isfinite(initial)));
    assert(dt > 0 && steps >= 1 && steps == floor(steps));
    if ~isempty(stance)
        assert(iscolumn(stance.waiting) && numel(stance.waiting) == numel(m));
        assert(all(isfinite(stance.waiting)) && all(stance.waiting >= 0 & stance.waiting <= 1));
        assert(numel(stance.parameters) == 7 && all(isfinite(stance.parameters)));
        assert(all(stance.parameters([1 5 7]) > 0) && all(stance.parameters(2:6) >= 0));
    end
    hasNeural = ~isempty(stance) && isfield(stance, 'neuralParameters');
    if hasNeural
        r = stance.neuralParameters;
        if numel(r) == 9
            r = [reshape(r, 1, 9), 0];
        end
        assert(numel(r) == 10 && all(isfinite(r)) && all(r(1:6) > 0));
        assert(r(7) >= 0 && r(7) <= 1 && all(r(8:10) >= 0));
        stance.neuralParameters = r;
    end
    divisive = hasNeural && isfield(stance, 'intentDecoder') && strcmp(stance.intentDecoder, 'divisive');
    softTonic = hasNeural && isfield(stance, 'intentDecoder') && strcmp(stance.intentDecoder, 'soft_tonic');
    if softTonic && useMex && exist('hill_neural_soft_tonic_mex', 'file') == 3
        [x, failed] = hill_neural_soft_tonic_mex(par, inputs, m, stance.waiting, initial, dt, steps, stance.parameters, stance.neuralParameters);
    elseif divisive && useMex && exist('hill_neural_divisive_mex', 'file') == 3
        [x, failed] = hill_neural_divisive_mex(par, inputs, m, stance.waiting, initial, dt, steps, stance.parameters, stance.neuralParameters);
    elseif hasNeural && ~divisive && ~softTonic && useMex && exist('hill_neural_mex', 'file') == 3
        [x, failed] = hill_neural_mex(par, inputs, m, stance.waiting, initial, dt, steps, stance.parameters, stance.neuralParameters);
    elseif ~hasNeural && ~isempty(stance) && useMex && exist('hill_sequence_mex', 'file') == 3
        [x, failed] = hill_sequence_mex(par, inputs, m, stance.waiting, initial, dt, steps, stance.parameters);
    elseif isempty(stance) && useMex && exist('hill_pair_mex', 'file') == 3
        [x, failed] = hill_pair_mex(par, inputs, m, initial, dt, steps);
    else
        [x, failed] = matlab_integrate(par, inputs, m, initial, dt, steps, stance);
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
    if ~isempty(stance)
        p.neuralCommands = x(:, 18:19);
        p.balanceWeight = x(:, 20);
        p.requestedHillTorque = x(:, 21);
        p.postureCorrectionRequest = x(:, 22);
        p.postureIntegral = x(:, 23);
        p.externalBalanceTorque = zeros(size(q));
        if hasNeural
            % Intent is an input-side neural estimate, not the predicted angle.
            p.motionIntent = x(:, 24);
        end
    end
end

function [out, failed] = matlab_integrate(par, a, m, initial, dt, steps, stance)
    front = [par(9), par(11), par(12), 0, par(25), par(10), par(26), par(27)];
    units = [par(1:8); front];
    q = initial(1);
    v = initial(2);
    h = dt / steps;
    n = size(a, 1);
    w = sin(par(28));
    gamma = reshape(units(:, 6), 1, 2);
    activationTau = reshape(units(:, 8), 1, 2);
    command = max(0, min(1, a(1, :)));
    blend = 0;
    integral = 0;
    desired = 0;
    correction = 0;
    control = [0 0];
    columns = 17;
    hasNeural = ~isempty(stance) && isfield(stance, 'neuralParameters');
    decoder = 'subtractive';
    if hasNeural && isfield(stance, 'intentDecoder'), decoder = stance.intentDecoder; end
    if ~isempty(stance)
        cfg = stance.parameters;
        blend = stance.waiting(1);
        columns = 23;
        if blend > 1e-8
            [control, desired, correction] = stance_command(par, units, cfg, q, v, m(1), integral);
        end
        if hasNeural
            r = reshape(stance.neuralParameters, 1, 10);
            fast = command(1);
            slow = command(1);
            neutral = excitation_targets(par, units, cfg(1), [max(0, -gravity(cfg(1), 0, par)), 0]);
            tonic = 0;
            if r(10) > 0
                loaded = excitation_targets(par, units, cfg(1), [max(0, -gravity(cfg(1), m(1), par)), 0]);
                tonic = max(0, loaded(1) - neutral(1));
            end
            signal = hill_motion_intent(command(1), tonic, r(10), decoder);
            intent = cfg(1) - r(1) * signal^r(2);
            columns = 24;
            if blend <= 1e-8
                [control, desired, correction] = movement_command(par, units, r, q, v, m(1), intent, command(2));
            end
            feedback = blend + (1 - blend) * r(7);
            command = (1 - feedback) * command + feedback * control;
        else
            command = (1 - blend) * command + blend * control;
        end
    end
    act = command.^gamma;
    fiberState = zeros(2, 6);
    oldl = zeros(1, 2);
    for u = 1:2
        fiberState(u, :) = fiber(units(u, :), q, act(u), 0, h, true, v, 3 - 2 * u, par(20 + u), par(22 + u), w);
        oldl(u) = sqrt(fiberState(u, 1)^2 + w^2);
    end
    out = zeros(n, columns);
    failed = false;
    for k = 1:n
        mom = 0;
        if k > 1
            for s = 1:steps
                f = s / steps;
                mass = m(k - 1) + f * (m(k) - m(k - 1));
                command = max(0, min(1, a(k - 1, :) + f * (a(k, :) - a(k - 1, :))));
                if ~isempty(stance)
                    targetBlend = max(0, min(1, stance.waiting(k - 1) + f * (stance.waiting(k) - stance.waiting(k - 1))));
                    blend = targetBlend + (blend - targetBlend) * exp(-h / cfg(7));
                    integral = max(-.5, min(.5, integral + h * blend * (cfg(1) - q)));
                    if hasNeural
                        fast = command(1) + (fast - command(1)) * exp(-h / .06);
                        slow = command(1) + (slow - command(1)) * exp(-h / .25);
                        tonic = 0;
                        if r(10) > 0
                            loaded = excitation_targets(par, units, cfg(1), [max(0, -gravity(cfg(1), mass, par)), 0]);
                            tonic = max(0, loaded(1) - neutral(1));
                        end
                        led = command(1) + r(9) * (fast - slow) / (.25 - .06);
                        intentSignal = hill_motion_intent(led, tonic, r(10), decoder);
                        intentTarget = cfg(1) - r(1) * intentSignal^r(2);
                        intent = intentTarget + (intent - intentTarget) * exp(-h / r(3));
                        post = [0 0];
                        moving = [0 0];
                        dp = 0; cp = 0; dm = 0; cm = 0;
                        if blend > 1e-8
                            [post, dp, cp] = stance_command(par, units, cfg, q, v, mass, integral);
                        end
                        if blend < 1 - 1e-8
                            [moving, dm, cm] = movement_command(par, units, r, q, v, mass, intent, command(2));
                        end
                        desired = blend * dp + (1 - blend) * dm;
                        correction = blend * cp + (1 - blend) * cm;
                        command = blend * post + (1 - blend) * ((1 - r(7)) * command + r(7) * moving);
                    elseif blend < 1e-8
                        blend = 0;
                        integral = integral * exp(-h);
                        desired = 0;
                        correction = 0;
                    else
                        [control, desired, correction] = stance_command(par, units, cfg, q, v, mass, integral);
                    end
                    if ~hasNeural
                        command = (1 - blend) * command + blend * control;
                    end
                end
                target = max(0, min(1, command)).^gamma;
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
        out(k, 1:17) = [q, v, force(1) - force(2), force', ce(1), pe(1), damp(1), ce(2), pe(2), damp(2), ...
                     oldl, act, max(abs(fiberState(:, 6))), mom];
        if ~isempty(stance)
            out(k, 18:23) = [command, blend, desired, correction, integral];
            if hasNeural
                out(k, 24) = intent;
            end
        end
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

function [command, desired, correction] = stance_command(par, units, cfg, q, v, mass, integral)
    % Convert a posture request into excitation; no extra balance torque acts.
    passiveTorque = zeros(1, 2);
    for u = 1:2
        fs = fiber(units(u, :), q, 0, 0, .001, true, 0, 3 - 2 * u, par(20 + u), par(22 + u), sin(par(28)));
        passiveTorque(u) = units(u, 1) * fs(2);
    end
    correction = cfg(5) * tanh((cfg(2) * (cfg(1) - q) - cfg(3) * v + cfg(4) * integral) / cfg(5));
    desired = -gravity(q, mass, par) + correction;
    difference = desired - passiveTorque(1) + passiveTorque(2);
    targets = passiveTorque + cfg(6) + [max(difference, 0), max(-difference, 0)];
    command = excitation_targets(par, units, q, targets);
end

function [command, desired, correction] = movement_command(par, units, r, q, v, mass, intent, frontInput)
    % Proprioceptive control changes excitation, never adds a joint torque.
    minimum = zeros(1, 2);
    for u = 1:2
        fs = fiber(units(u, :), q, 0, 0, .001, true, 0, 3 - 2 * u, par(20 + u), par(22 + u), sin(par(28)));
        minimum(u) = units(u, 1) * fs(2);
    end
    fs = fiber(units(2, :), q, max(0, min(1, frontInput))^units(2, 6), 0, .001, true, 0, -1, par(22), par(24), sin(par(28)));
    baseFront = units(2, 1) * fs(2);
    correction = r(4) * (intent - q) - r(5) * v;
    request = -gravity(q, mass, par) + correction;
    desired = request;
    if request < 0
        desired = r(6) * tanh(request / r(6));
    end
    targets = [max(minimum(1), desired + baseFront), max(baseFront, minimum(1) - desired)] + r(8);
    command = excitation_targets(par, units, q, targets);
end

function command = excitation_targets(par, units, q, targets)
    % Static inversion requests excitation; dynamic Hill forces remain free.
    command = zeros(1, 2);
    for u = 1:2
        unit = units(u, :);
        lo = par(20 + u) / unit(2);
        w = sin(par(28));
        slack = par(22 + u) - lo * sqrt(unit(3)^2 - w^2);
        lmt = par(22 + u) - (3 - 2 * u) * par(20 + u) * (q - 1.4);
        ft = max(0, targets(u) / unit(1));
        lt = slack * (1 + unit(7) * sqrt(ft));
        x = (lmt - lt) / lo;
        if x <= 0
            continue
        end
        l = sqrt(x^2 + w^2);
        activation = (ft / (x / l) - passive(l)) / max(exp(-((l - 1) / .45)^2), 1e-8);
        command(u) = max(0, min(1, activation))^(1 / unit(6));
    end
end
