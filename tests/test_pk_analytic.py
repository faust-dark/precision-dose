"""Analytic correctness of the compartment engine.

Two independent oracles are used:
- the closed-form biexponential solution for a two-compartment IV bolus,
- the closed-form steady-state solution for one-compartment oral dosing.

If the matrix exponential, the infusion superposition or the cycle fixed point
is wrong, these tests fail with a numerical difference rather than a vague
"output changed".
"""

import math
import unittest

from precision_dose import pk


class TestTwoCompartment(unittest.TestCase):
    def setUp(self):
        self.cl, self.v1, self.v2, self.q = 4.4, 37.9, 44.2, 6.95

    def analytic(self, dose, t):
        k10 = self.cl / self.v1
        k12 = self.q / self.v1
        k21 = self.q / self.v2
        s = k10 + k12 + k21
        disc = math.sqrt(s * s - 4 * k10 * k21)
        alpha, beta = (s + disc) / 2, (s - disc) / 2
        a = dose / self.v1 * (alpha - k21) / (alpha - beta)
        b = dose / self.v1 * (k21 - beta) / (alpha - beta)
        return a * math.exp(-alpha * t) + b * math.exp(-beta * t)

    def test_bolus_matches_biexponential(self):
        model = pk.LinearPK(self.cl, self.v1, self.v2, self.q)
        for t in (0.5, 2.0, 8.0, 24.0, 48.0):
            got = model.simulate([(0.0, 1000.0, 0.0001, "iv")], [t])[0][1]
            self.assertAlmostEqual(got, self.analytic(1000.0, t), delta=1e-6 * max(1, got))

    def test_half_life_is_terminal_phase(self):
        model = pk.LinearPK(self.cl, self.v1, self.v2, self.q)
        thalf = model.half_life()
        # sample well past the distribution phase (alpha t1/2 is under an hour here)
        c0 = self.analytic(1000.0, 100.0)
        c1 = self.analytic(1000.0, 100.0 + thalf)
        self.assertAlmostEqual(c1 / c0, 0.5, 4)


class TestMassBalance(unittest.TestCase):
    def test_iv_steady_state_auc_equals_dose_over_cl(self):
        model = pk.LinearPK(4.0, 35.0, 50.0, 6.0)
        for tau in (8.0, 12.0, 24.0, 48.0):
            ev = model.evaluate(pk.Regimen([(0.0, 1000.0)], tau, 1.0, "iv"), grid_h=0.2)
            daily = 1000.0 * 24.0 / tau
            self.assertAlmostEqual(ev["auc24"], daily / 4.0, delta=1e-6 * daily)

    def test_oral_auc_equals_f_times_dose_over_cl(self):
        model = pk.LinearPK(1.1, 45.0, ka=1.2, f=0.95, lag_h=0.5)
        ev = model.evaluate(pk.Regimen([(0.0, 10.0)], 24.0, 0.0, "oral"), grid_h=0.05)
        self.assertAlmostEqual(ev["auc24"], 0.95 * 10.0 / 1.1, delta=1e-6)

    def test_amounts_never_negative_and_bounded(self):
        model = pk.LinearPK(2.0, 30.0, 60.0, 5.0)
        ev = model.evaluate(pk.Regimen([(0.0, 750.0)], 12.0, 2.0, "iv"), grid_h=0.1)
        amounts = [c for _, c in ev["curve"]]
        self.assertGreaterEqual(min(amounts), 0.0)
        self.assertLess(max(amounts), 750.0 / 30.0 * 2.0)

    def test_steady_state_is_a_fixed_point(self):
        model = pk.LinearPK(3.0, 40.0, ka=2.0, f=1.0)
        reg = pk.Regimen([(0.0, 500.0)], 12.0, 0.0, "oral")
        x0, segs = model.steady_state(reg)
        x_end = list(x0)
        for h, rate, extra, bolus in segs:
            for idx, amt in bolus.items():
                x_end[idx] += amt
            x_end = model._advance(x_end, rate, h, extra)
        for a, b in zip(x0, x_end):
            self.assertAlmostEqual(a, b, delta=1e-9 * max(1.0, abs(a)))


class TestOneCompartment(unittest.TestCase):
    def test_half_life(self):
        model = pk.LinearPK(2.0, 40.0)
        self.assertAlmostEqual(model.half_life(), math.log(2) * 40.0 / 2.0, 9)

    def test_trough_closed_form(self):
        cl, v, tau, dose = 2.0, 40.0, 12.0, 500.0
        model = pk.LinearPK(cl, v)
        ev = model.evaluate(pk.Regimen([(0.0, dose)], tau, 0.0, "iv"), grid_h=0.05)
        k = cl / v
        expected = -math.log(dose / v * math.exp(-k * tau) / (1 - math.exp(-k * tau))) * 0
        expected = (dose / v) * math.exp(-k * tau) / (1 - math.exp(-k * tau))
        self.assertAlmostEqual(ev["trough"], expected, delta=1e-6 * expected)


class TestMichaelisMenten(unittest.TestCase):
    def test_steady_state_relation(self):
        model = pk.MichaelisMentenPK(415.0, 5.7, 56.0)
        css = model.steady_state_concentration(300.0)
        # R = Vmax*C/(Km+C) must hold exactly at the returned concentration
        self.assertAlmostEqual(300.0, 415.0 * css / (5.7 + css), 9)

    def test_no_steady_state_above_vmax(self):
        model = pk.MichaelisMentenPK(415.0, 5.7, 56.0)
        self.assertIsNone(model.steady_state_concentration(500.0))

    def test_half_saturation_gives_km(self):
        model = pk.MichaelisMentenPK(400.0, 5.0, 56.0)
        self.assertAlmostEqual(model.steady_state_concentration(200.0), 5.0, 6)


if __name__ == "__main__":
    unittest.main()
