"""Compartmental PK models.

Linear models are solved exactly: the state transition over a step h is the
matrix exponential of the rate matrix, and a periodic regimen is evaluated at
its fixed point (one linear solve) rather than by iterating to convergence.
Mass balance therefore holds to machine precision and there is no step-size
error to defend.

The only model integrated numerically is phenytoin's Michaelis-Menten
elimination, where no closed form exists; RK4 with a reported step is used for
the curve, while the steady-state target is evaluated from the exact algebraic
relation R = Vmax*Css/(Km+Css).
"""

import math

from . import matmath
from .errors import E_BAD_VALUE, EngineError

DEFAULT_GRID_H = 0.25


class Regimen:
    """One cycle of a dosing regimen.

    doses: list of (time_within_cycle_h, amount_in_drug_units)
    """

    def __init__(self, doses, cycle_h, infusion_h=0.0, route="iv"):
        self.doses = sorted((float(t), float(a)) for t, a in doses)
        self.cycle_h = float(cycle_h)
        self.infusion_h = float(infusion_h)
        self.route = route

    def as_dict(self):
        return {
            "doses": [
                {"time_h": round(t, 4), "amount": round(a, 6)} for t, a in self.doses
            ],
            "cycle_h": round(self.cycle_h, 4),
            "infusion_h": round(self.infusion_h, 4),
            "route": self.route,
        }


class LinearPK:
    """One- or two-compartment linear model with optional first-order absorption.

    Volumes and clearances are *apparent* when the entry declares oral dosing:
    F is folded into the parameters, which is the standard popPK convention for
    drugs monitored after oral administration (tacrolimus).
    """

    def __init__(self, cl, v1, v2=None, q=None, ka=None, f=1.0, lag_h=0.0, conc_factor=1.0):
        self.conc_factor = float(conc_factor)
        self.cl = float(cl)
        self.v1 = float(v1)
        self.v2 = float(v2) if v2 else None
        self.q = float(q) if q else None
        self.ka = float(ka) if ka else None
        self.f = float(f)
        self.lag_h = float(lag_h)
        if self.cl <= 0 or self.v1 <= 0:
            raise EngineError(E_BAD_VALUE, "CL and V1 must be positive")
        self.oral = self.ka is not None
        self.n = 1 + (1 if self.v2 else 0) + (1 if self.oral else 0)
        self.central = 1 if self.oral else 0
        if self.oral:
            self.depot = 0

    # -- matrices ---------------------------------------------------------
    def matrix(self, cl_extra=0.0):
        cl = self.cl + cl_extra
        n = self.n
        m = [[0.0] * n for _ in range(n)]
        c = self.central
        if self.oral:
            m[0][0] = -self.ka
            # the dose entering the depot is already multiplied by F; applying F
            # again here would double-count the bioavailability
            m[c][0] = self.ka
        if self.v2:
            m[c][c] = -(cl + self.q) / self.v1
            m[c][c + 1] = self.q / self.v2
            m[c + 1][c] = self.q / self.v1
            m[c + 1][c + 1] = -self.q / self.v2
        else:
            m[c][c] = -cl / self.v1
        return m

    # -- exact step -------------------------------------------------------
    def _advance(self, x, rate, h, cl_extra):
        if h <= 0:
            return list(x)
        m = self.matrix(cl_extra)
        e = matmath.mat_exp(matmath.mat_scale(m, h))
        x = matmath.mat_vec(e, x)
        if rate:
            u = [0.0] * self.n
            u[self.central] = rate
            delta = [a - b for a, b in zip(matmath.mat_vec(e, u), u)]
            try:
                phi = matmath.solve(m, delta)
            except ZeroDivisionError:
                phi = [ui * h for ui in u]
            x = [a + b for a, b in zip(x, phi)]
        return x

    # -- timeline ---------------------------------------------------------
    def timeline(self, regimen, dialysis=None, hd_cl_h=0.0):
        """Piecewise-constant segments covering exactly one cycle.

        Each segment: (duration_h, rate_into_central, extra_clearance, bolus).
        bolus is a dict {state_index: amount} applied at the segment start.
        """
        points = {0.0, regimen.cycle_h}
        infusions = []
        boluses = []
        hd_sessions = []

        for t, amt in regimen.doses:
            if regimen.route == "oral":
                t_bolus = t + self.lag_h
                boluses.append((t_bolus, amt * self.f))
                points.add(t_bolus)
            elif regimen.infusion_h <= 0:
                # a zero-duration infusion is a bolus into the central compartment
                boluses.append((t, amt))
                points.add(t)
            else:
                start, end = t, t + regimen.infusion_h
                rate = amt / regimen.infusion_h
                infusions.append((start, end, rate))
                points.update((start, end))

        if dialysis and dialysis.get("modality") == "hd":
            dur = float(dialysis["hours_per_session"])
            for i in range(int(dialysis["sessions_per_week"])):
                start = float(dialysis["first_session_h"]) + i * float(
                    dialysis["day_interval_h"]
                )
                hd_sessions.append((start, start + dur))
                points.update((start, start + dur))

        # clip to the cycle
        points = sorted(p for p in points if 0.0 <= p <= regimen.cycle_h + 1e-12)
        if points[0] > 0.0:
            points.insert(0, 0.0)
        if points[-1] < regimen.cycle_h:
            points.append(regimen.cycle_h)

        segs = []
        for a, b in zip(points, points[1:]):
            h = b - a
            if h <= 1e-12:
                continue
            mid = 0.5 * (a + b)
            rate = sum(r for s, e, r in infusions if s <= mid < e)
            extra = 0.0
            if any(s <= mid < e for s, e in hd_sessions):
                extra = hd_cl_h
            bolus = {}
            for t_b, amt in boluses:
                if abs(t_b - a) < 1e-9:
                    idx = self.depot if self.oral else self.central
                    bolus[idx] = bolus.get(idx, 0.0) + amt
            segs.append((h, rate, extra, bolus))
        return segs

    def _compose(self, segs):
        n = self.n
        a = matmath.identity(n)
        b = [0.0] * n
        for h, rate, cl_extra, bolus in segs:
            # a bolus lands at the START of its segment, so it must be inside
            # the transition: x_out = E(x_in + B) + phi, not E x_in + phi + B
            if bolus:
                for idx, amt in bolus.items():
                    b[idx] += amt
            m = self.matrix(cl_extra)
            e = matmath.mat_exp(matmath.mat_scale(m, h))
            b = matmath.mat_vec(e, b)
            if rate:
                u = [0.0] * n
                u[self.central] = rate
                delta = [x - y for x, y in zip(matmath.mat_vec(e, u), u)]
                try:
                    phi = matmath.solve(m, delta)
                except ZeroDivisionError:
                    phi = [ui * h for ui in u]
                b = [x + y for x, y in zip(b, phi)]
            a = matmath.mat_mul(e, a)
        return a, b

    def steady_state(self, regimen, dialysis=None, hd_cl_h=0.0):
        segs = self.timeline(regimen, dialysis, hd_cl_h)
        a, b = self._compose(segs)
        n = self.n
        i_minus_a = [
            [(1.0 if i == j else 0.0) - a[i][j] for j in range(n)] for i in range(n)
        ]
        x0 = matmath.solve(i_minus_a, b)
        return x0, segs

    # -- evaluation --------------------------------------------------------
    def evaluate(
        self, regimen, dialysis=None, hd_cl_h=0.0, grid_h=DEFAULT_GRID_H
    ):
        x0, segs = self.steady_state(regimen, dialysis, hd_cl_h)
        x = list(x0)
        t = 0.0
        curve = [(0.0, self.conc_factor * x[self.central] / self.v1)]
        auc = 0.0
        for h, rate, cl_extra, bolus in segs:
            if bolus:
                for idx, amt in bolus.items():
                    x[idx] += amt
                curve.append((t, self.conc_factor * x[self.central] / self.v1))
            m = self.matrix(cl_extra)
            steps = max(1, int(math.ceil(h / grid_h)))
            dt = h / steps
            u = [0.0] * self.n
            if rate:
                u[self.central] = rate
            for _ in range(steps):
                x_next = self._advance(x, rate, dt, cl_extra)
                delta = [x_next[i] - x[i] - u[i] * dt for i in range(self.n)]
                try:
                    integral = matmath.solve(m, delta)
                except ZeroDivisionError:
                    integral = [d * dt / 2.0 for d in delta]
                auc += self.conc_factor * integral[self.central] / self.v1
                x = x_next
                t += dt
                curve.append((t, self.conc_factor * x[self.central] / self.v1))
        cmax = max(c for _, c in curve)
        cmin = min(c for _, c in curve)
        tmax = max(curve, key=lambda p: (p[1], -p[0]))[0]
        # mean daily exposure: exact for any cycle length (regimen.cycle_h)
        auc24 = auc * (24.0 / regimen.cycle_h) if regimen.cycle_h else auc
        return {
            "curve": curve,
            "auc_cycle": auc,
            "auc24": auc24,
            "cmax": cmax,
            "cmin": cmin,
            "tmax_h": tmax,
            "trough": curve[-1][1],
            "state_at_cycle_start": x0,
            "cycle_h": regimen.cycle_h,
            "grid_h": grid_h,
        }

    def effective_half_life(self, regimen, points=200):
        """Accumulation-based effective half-life over the dosing interval.

        For a multi-compartment drug fitted to trough data the terminal phase is
        often unidentifiable (the peripheral volume is a sink, not a measurement),
        so the terminal half-life can be a number nobody should act on. The
        effective half-life is defined from what the model does predict:
        tau * ln2 / ln(AUC_ss / AUC_single) over one interval.
        """
        if regimen.cycle_h <= 0:
            return None
        times = [i * regimen.cycle_h / points for i in range(points + 1)]
        doses = [
            (t, a, regimen.infusion_h if regimen.route != "oral" else 0.0, regimen.route)
            for t, a in regimen.doses
        ]
        single = self.simulate(doses, times)
        auc_single = 0.0
        for (t0, c0), (t1, c1) in zip(single, single[1:]):
            auc_single += 0.5 * (c0 + c1) * (t1 - t0)
        ev = self.evaluate(regimen)
        auc_ss = ev["auc_cycle"]
        if auc_single <= 0 or auc_ss <= auc_single:
            return None
        # accumulation ratio -> half-life of the one-compartment model that would
        # accumulate identically: t = tau*ln2 / ln(R/(R-1))
        ratio = auc_ss / auc_single
        if ratio <= 1.0 + 1e-12:
            return None
        return regimen.cycle_h * math.log(2.0) / math.log(ratio / (ratio - 1.0))

    def half_life(self):
        if self.v2:
            k10 = self.cl / self.v1
            k12 = self.q / self.v1
            k21 = self.q / self.v2
            s = k10 + k12 + k21
            disc = max(s * s - 4.0 * k10 * k21, 0.0)
            beta = (s - math.sqrt(disc)) / 2.0
            return math.log(2.0) / beta
        return math.log(2.0) * self.v1 / self.cl

    # -- targeted sampling (used by PTA, where a full grid is wasteful) -----
    def sample(self, regimen, times, dialysis=None, hd_cl_h=0.0):
        """Concentration at specific times within the cycle, at steady state."""
        x0, segs = self.steady_state(regimen, dialysis, hd_cl_h)
        want = sorted(float(t) for t in times)
        out = []
        x = list(x0)
        t = 0.0
        idx = 0
        for h, rate, extra, bolus in segs:
            if bolus:
                for i, amt in bolus.items():
                    x[i] += amt
            while idx < len(want) and want[idx] <= t + 1e-9:
                out.append((want[idx], self.conc_factor * x[self.central] / self.v1))
                idx += 1
            seg_end = t + h
            while idx < len(want) and want[idx] < seg_end - 1e-12:
                target = want[idx]
                x = self._advance(x, rate, target - t, extra)
                t = target
                out.append((t, self.conc_factor * x[self.central] / self.v1))
                idx += 1
            x = self._advance(x, rate, seg_end - t, extra)
            t = seg_end
        while idx < len(want):
            out.append((want[idx], self.conc_factor * x[self.central] / self.v1))
            idx += 1
        return out

    def simulate(self, doses, times):
        """Concentration after an arbitrary dosing history starting from zero.

        doses: list of (start_h, amount, infusion_min, route)
        times: sampling times relative to the same origin
        Used by the Bayesian estimator, where the history is not periodic.
        """
        points = {0.0}
        infusions = []
        boluses = []
        for start, amount, infusion_min, route in doses:
            start = float(start)
            if route == "oral":
                t_b = start + self.lag_h
                boluses.append((t_b, float(amount) * self.f))
                points.add(t_b)
            else:
                end = start + float(infusion_min) / 60.0
                infusions.append((start, end, float(amount) / max(end - start, 1e-9)))
                points.update((start, end))
        horizon = max(list(points) + [float(t) for t in times]) + 1e-6
        points = sorted(p for p in points if 0.0 <= p <= horizon)
        if points[0] > 0:
            points.insert(0, 0.0)
        if points[-1] < horizon:
            points.append(horizon)

        x = [0.0] * self.n
        t = 0.0
        out = []
        want = sorted(float(s) for s in times)
        idx = 0
        for a, b in zip(points, points[1:]):
            for t_b, amt in boluses:
                if abs(t_b - a) < 1e-9:
                    idxb = self.depot if self.oral else self.central
                    x[idxb] += amt
            rate = sum(r for s, e, r in infusions if s <= a + 0.5 * (b - a) < e)
            while idx < len(want) and want[idx] <= b + 1e-12:
                target = max(want[idx], t)
                x = self._advance(x, rate, target - t, 0.0)
                t = target
                out.append((want[idx], self.conc_factor * x[self.central] / self.v1))
                idx += 1
            x = self._advance(x, rate, b - t, 0.0)
            t = b
        return out



class MichaelisMentenPK:
    """Capacity-limited elimination (phenytoin): R = Vmax*C/(Km+C)."""

    def __init__(self, vmax_mg_day, km_mg_L, vd_L, f=1.0, conc_factor=1.0):
        self.conc_factor = float(conc_factor)
        self.vmax = float(vmax_mg_day) / 24.0
        self.vmax_day = float(vmax_mg_day)
        self.km = float(km_mg_L)
        self.vd = float(vd_L)
        self.f = float(f)
        if self.vmax <= 0 or self.km <= 0 or self.vd <= 0:
            raise EngineError(E_BAD_VALUE, "Vmax, Km and Vd must be positive")

    def steady_state_concentration(self, dose_per_day, f=None):
        """Exact algebraic steady state; None when input exceeds capacity."""
        f = self.f if f is None else f
        r = f * dose_per_day
        if r >= self.vmax_day:
            return None
        return self.km * r / (self.vmax_day - r)

    def evaluate(self, regimen, grid_h=1.0, horizon_h=720.0):
        doses = regimen.doses
        cycle = regimen.cycle_h
        steps_per_h = max(1, int(round(1.0 / grid_h)))
        dt = 1.0 / steps_per_h
        horizon = max(horizon_h, 25 * cycle)

        def deriv(a):
            c = max(a / self.vd, 0.0)
            return -(self.vmax * c / (self.km + c)) if c > 0 else 0.0

        a = 0.0
        t = 0.0
        curve = [(0.0, 0.0)]
        n_steps = int(horizon * steps_per_h)
        for k in range(n_steps):
            phase = t % cycle
            for td, amt in doses:
                if abs(phase - td) < dt * 0.5 and t > 0:
                    a += self.f * amt
            k1 = deriv(a)
            k2 = deriv(a + 0.5 * dt * k1)
            k3 = deriv(a + 0.5 * dt * k2)
            k4 = deriv(a + dt * k3)
            a = a + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            t = (k + 1) * dt
            curve.append((t, self.conc_factor * a / self.vd))
        tail = [p for p in curve if p[0] >= horizon - cycle]
        cmax = max(c for _, c in tail)
        cmin = min(c for _, c in tail)
        auc = 0.0
        for (t0, c0), (t1, c1) in zip(tail, tail[1:]):
            auc += 0.5 * (c0 + c1) * (t1 - t0)
        return {
            "curve": curve,
            "auc_cycle": auc,
            "auc24": auc * (24.0 / cycle) if cycle else auc,
            "cmax": cmax,
            "cmin": cmin,
            "tmax_h": max(tail, key=lambda p: (p[1], -p[0]))[0],
            "trough": tail[-1][1],
            "css_avg": auc / cycle if cycle else None,
            "integration_step_h": dt,
        }
