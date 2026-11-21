"""Mathematical CEINMS stiff-tendon force audit, with explicit component forces.

Equations and spline boundary conditions were checked against CEINMS-RT core
commit 5e53f3293ff4a5662ce575e83d2610c146de35ff. This is an independent replay
of exported states, not the acquisition executable or a ground-truth calibration.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np
from scipy.interpolate import PPoly


class MuscleCurve:
    """Cubic interpolant with endpoint third derivatives from divided differences."""

    def __init__(self, x, y):
        self.x, self.y = np.asarray(x), np.asarray(y)
        h = np.diff(x)
        slope = np.diff(y) / h
        n = len(x)
        system = np.zeros((n, n))
        rhs = np.zeros(n)
        system[0, :2] = [-h[0], h[0]]
        system[-1, -2:] = [h[-1], -h[-1]]
        for k in range(1, n - 1):
            system[k, k - 1 : k + 2] = [h[k - 1], 2 * (h[k - 1] + h[k]), h[k]]
            rhs[k] = slope[k] - slope[k - 1]
        if n > 3:
            rhs[0] = (
                (
                    (slope[2] - slope[1]) / (x[3] - x[1])
                    - (slope[1] - slope[0]) / (x[2] - x[0])
                )
                * h[0] ** 2
                / (x[3] - x[0])
            )
            rhs[-1] = (
                -(
                    (slope[-1] - slope[-2]) / (x[-1] - x[-3])
                    - (slope[-2] - slope[-3]) / (x[-2] - x[-4])
                )
                * h[-1] ** 2
                / (x[-1] - x[-4])
            )
        coefficients = np.linalg.solve(system, rhs)
        linear = slope - h * (2 * coefficients[:-1] + coefficients[1:])
        cubic = np.diff(coefficients) / h
        self.spline = PPoly(np.array([cubic, 3 * coefficients[:-1], linear, y[:-1]]), x)

    def __call__(self, values):
        return self.spline(np.clip(values, self.x[0], self.x[-1]))

    def derivative(self, values):
        inside = (values >= self.x[0]) & (values <= self.x[-1])
        return np.where(inside, self.spline.derivative()(values), 0)


class StiffMuscles:
    """Use the subject's saved parameters in its exported muscle column order."""

    def __init__(self, xml, names):
        tree = ET.parse(xml)
        muscles = {m.findtext("name"): m for m in tree.iter("muscle")}
        defaults = tree.find("muscleDefault")
        self.names = tuple(names)

        def parameters(key):
            return np.array([float(muscles[n].findtext(key)) for n in names])

        self.l0 = parameters("optimalFiberLength")
        self.slack = parameters("tendonSlackLength")
        self.pennation_parameter = parameters("pennationAngle")
        self.fmax = parameters("maxIsometricForce") * parameters("strengthCoefficient")
        self.old_shape = parameters("shapeFactor")
        self.length_change = float(defaults.findtext("percentageChange"))
        self.damping = float(defaults.findtext("damping"))
        self.curves = {}
        for node in defaults.findall("Curve"):
            x, y = [
                np.fromstring(node.findtext(k), sep=" ") for k in ("xPoints", "yPoints")
            ]
            self.curves[node.findtext("name")] = MuscleCurve(x, y)
        self.groups = {
            channel.findtext("name"): [
                self.names.index(n)
                for n in channel.findtext("muscleSequence").split()
                if n in self.names
            ]
            for channel in tree.iter("Channel")
        }

    def lengths(self, lmt, activation):
        # The saved runtime logs match DEGREE interpretation of this parameter;
        # this is a replay convention, not a claim about its anatomical units.
        effective_l0 = self.l0 * (1 + self.length_change * (1 - activation))
        thickness = effective_l0 * np.sin(np.deg2rad(self.pennation_parameter))
        return np.sqrt(thickness**2 + (lmt - self.slack) ** 2)

    def forces(self, activation, length, velocity, pennation_degrees=None):
        normalized_length = length / self.l0
        adjusted_length = length / (
            self.l0 * (1 + self.length_change * (1 - activation))
        )
        normalized_velocity = np.clip(velocity / self.l0, -10, 10) / 10
        if pennation_degrees is None:
            sine = self.l0 * np.sin(np.deg2rad(self.pennation_parameter)) / length
            projection = np.sqrt(1 - np.clip(sine, 0, 1) ** 2)
        else:
            projection = np.cos(np.deg2rad(pennation_degrees))
        multiplier = self.fmax * projection
        active = (
            multiplier
            * activation
            * self.curves["activeForceLength"](adjusted_length)
            * self.curves["forceVelocity"](normalized_velocity)
        )
        elastic = multiplier * self.curves["passiveForceLength"](normalized_length)
        damping = multiplier * self.damping * normalized_velocity
        return dict(
            active=active,
            passive=elastic + damping,
            elastic=elastic,
            damping=damping,
            total=np.maximum(0, active + elastic + damping),
            activeOnly=np.maximum(0, active),
        )

    def recover_recorded_activation(self, lmt, length, recorded_active, initial):
        """Recover the asynchronous state, not physiological parameters.

        Back-channel activation is determined from the recorded stiff-tendon
        thickness. Zero-pennation abdominal channels use their active-force
        equation. The result is used only for a same-state counterfactual; it
        is not an independent forward validation against measured force.
        """
        activation = np.zeros_like(length)
        denominator = self.l0 * np.sin(np.deg2rad(self.pennation_parameter))
        thickness = np.sqrt(
            np.maximum(
                (length - abs(lmt - self.slack)) * (length + abs(lmt - self.slack)), 0
            )
        )
        for _, indices in self.groups.items():
            if not indices:
                continue
            supported = [k for k in indices if denominator[k] > 1e-6]
            if supported:
                recovered = (
                    1
                    + self.length_change
                    - thickness[:, supported] / denominator[supported]
                ) / self.length_change
                value = np.median(recovered, axis=1)
            else:
                # Pick the member whose active curve carries the most force.
                k = max(indices, key=lambda j: np.mean(abs(recorded_active[:, j])))
                value = initial[:, k].copy()
                force = recorded_active[:, k] / self.fmax[k]
                for _ in range(12):
                    factor = 1 + self.length_change * (1 - value)
                    q = length[:, k] / (self.l0[k] * factor)
                    fl = self.curves["activeForceLength"](q)
                    gradient = (
                        fl
                        + value
                        * self.curves["activeForceLength"].derivative(q)
                        * q
                        * self.length_change
                        / factor
                    )
                    change = np.divide(
                        value * fl - force,
                        gradient,
                        out=np.zeros_like(value),
                        where=abs(gradient) > 1e-9,
                    )
                    value -= change
                if not np.isfinite(value).all() or np.max(abs(value)) > 2:
                    raise ValueError("Activation-state inversion failed")
            activation[:, indices] = value[:, None]
        return activation


def causal_velocity(times, length, history=15):
    """Endpoint derivative of a rolling cubic, matching the online trace rule."""
    velocity = np.zeros_like(length)
    times = times - times[0]
    for end in range(1, len(times)):
        lo = max(0, end - history + 1)
        x, y = times[lo : end + 1], length[lo : end + 1]
        h = np.diff(x)
        slope = np.diff(y, axis=0) / h[:, None]
        n = len(x)
        if n == 2:
            velocity[end] = slope[0]
            continue
        system = np.zeros((n, n))
        rhs = np.zeros((n, y.shape[1]))
        system[0, :2] = [-h[0], h[0]]
        system[-1, -2:] = [h[-1], -h[-1]]
        for k in range(1, n - 1):
            system[k, k - 1 : k + 2] = [h[k - 1], 2 * (h[k - 1] + h[k]), h[k]]
            rhs[k] = slope[k] - slope[k - 1]
        if n > 3:
            rhs[0] = (
                (
                    (slope[2] - slope[1]) / (x[3] - x[1])
                    - (slope[1] - slope[0]) / (x[2] - x[0])
                )
                * h[0] ** 2
                / (x[3] - x[0])
            )
            rhs[-1] = (
                -(
                    (slope[-1] - slope[-2]) / (x[-1] - x[-3])
                    - (slope[-2] - slope[-3]) / (x[-2] - x[-4])
                )
                * h[-1] ** 2
                / (x[-1] - x[-4])
            )
        coefficients = np.linalg.solve(system, rhs)
        velocity[end] = slope[-1] + h[-1] * (coefficients[-2] + 2 * coefficients[-1])
    return velocity
