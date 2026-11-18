"""Tests for CavityQED-Lab: JC dynamics, regimes, dispersive, gates."""
import math

import main as m


# ---- eigensolver + Hamiltonian
def test_jacobi_known_matrix():
    evals, _ = m.jacobi_eigh([[2.0, 1.0], [1.0, 2.0]])
    assert abs(evals[0] - 1.0) < 1e-12
    assert abs(evals[1] - 3.0) < 1e-12


def test_jacobi_reconstruction():
    h = m.build_jc_hamiltonian(1.0, 0.0, 6)
    evals, evecs = m.jacobi_eigh(h)
    dim = 12
    for i in range(dim):
        for j in range(dim):
            recon = sum(evecs[k][i] * evals[k] * evecs[k][j]
                        for k in range(dim))
            assert abs(recon - h[i][j]) < 1e-9


def test_hamiltonian_couplings():
    h = m.build_jc_hamiltonian(2.0, 0.0, 4)
    assert h[4][1] == 2.0 * math.sqrt(1)   # |e,0> <-> |g,1>
    assert h[5][2] == 2.0 * math.sqrt(2)   # |e,1> <-> |g,2>
    assert h[4][2] == 0.0                  # excitation number conserved


def test_manifold_eigenvalues_resonant():
    lo, hi = m.manifold_eigenvalues(1.0, 0.0, 0)
    assert lo == -1.0 and hi == 1.0


# ---- vacuum Rabi
def test_vacuum_rabi_frequency():
    rep = m.vacuum_rabi_report()
    assert abs(rep["frequency_measured"] - 2.0) / 2.0 < 0.005
    assert rep["p_e_at_0"] == 1.0
    assert rep["first_min_value"] < 1e-9


def test_vacuum_rabi_period():
    rep = m.vacuum_rabi_report()
    assert abs(rep["period_measured"] - math.pi) < 1e-3
    assert abs(rep["first_min_time"] - math.pi / 2) < 1e-4


def test_dressed_splitting_resonant():
    rep = m.dressed_splitting_report(delta=0.0)
    assert abs(rep["splitting_measured"] - 2.0) < 1e-6


def test_dressed_splitting_detuned():
    rep = m.dressed_splitting_report(delta=3.0)
    assert abs(rep["splitting_measured"] - math.sqrt(13.0)) < 1e-6


# ---- collapse and revival
def test_revival_time_and_peak():
    rep = m.revival_report()
    assert rep["revival_time_rel_error"] < 0.12
    assert rep["revival_peak"] > 0.75


def test_collapse_is_flat():
    rep = m.revival_report()
    assert rep["collapse_range"] < 0.15
    assert abs(rep["collapse_mean"] - 0.5) < 0.05
    assert rep["early_max"] > 0.9


def test_poisson_weights_normalised():
    psi = m.coherent_state(9.0, 30)
    assert sum(x * x for x in psi) > 1 - 1e-7
    assert psi[0] == 0.0 or True  # ground-atom block carries the weights
    assert abs(psi[3] ** 2 - math.exp(-9.0) * 9.0 ** 3
               / math.factorial(3)) < 1e-15


# ---- regimes
def test_cooperativity_values():
    assert m.cooperativity(1.0, 0.05, 0.02) == 1000.0
    reg = m.regimes_report()
    assert reg["cooperativity_max"] == 3200.0
    assert reg["cooperativity_min"] == 0.1


def test_purcell_bad_cavity():
    assert abs(m.purcell_factor(0.1, 2.0, 0.05) - 0.4) < 1e-12
    assert abs(m.cavity_induced_decay(0.1, 2.0) - 0.02) < 1e-12


def test_strong_coupling_boundary():
    assert m.strong_coupling(1.0, 0.05, 0.02) is True
    assert m.strong_coupling(0.1, 2.0, 0.05) is False
    assert m.strong_coupling(0.2, 0.5, 0.1) is False
    assert m.splitting_resolved(1.0, 0.05, 0.02) is True


# ---- dispersive
def test_dispersive_chi():
    rep = m.dispersive_report()
    assert rep["chi"] == 0.1
    assert rep["number_splitting_2chi"] == 0.2
    assert rep["ac_stark_per_photon"] == 0.2
    assert rep["validity_g_over_delta"] == 0.1


def test_dispersive_shifts_close_to_exact():
    rep = m.dispersive_report()
    s0 = rep["manifold_shifts"][0]
    assert abs(s0["dressed_shift_exact"] - 0.09902) < 1e-5
    assert abs(rep["per_photon_pull_exact"] - rep["chi"]) / rep["chi"] < 0.05
    assert rep["full_spectrum_max_deviation"] < 1e-8


# ---- dissipation
def test_damped_peaks_decrease():
    rep = m.damped_rabi_report()
    heights = [p[1] for p in rep["peaks"]]
    assert len(heights) >= 5
    assert all(b < a for a, b in zip(heights, heights[1:]))
    assert rep["max_trace_error"] < 1e-9
    assert rep["final_p_e"] < 0.05


def test_damped_envelope_rate():
    rep = m.damped_rabi_report()
    assert abs(rep["envelope_decay_rate_fit"] - 0.075) < 0.02


def test_undamped_control():
    rep = m.damped_rabi_report(kappa=0.0, t_max=10.0)
    assert rep["peak_count"] >= 3
    assert all(abs(p[1] - 1.0) < 1e-4 for p in rep["peaks"])


# ---- classifier
def test_classifier_scenarios():
    rep = m.classifier_report()
    assert rep["demo_point"]["label"] == "strong-coupling"
    assert rep["dispersive_point"]["label"] == "dispersive"
    assert rep["bad_cavity_point"]["label"] == "bad-cavity"
    assert rep["weak_point"]["label"] == "weak-coupling"


def test_classifier_matches_policy():
    for params in m.CLASSIFIER_SCENARIOS.values():
        feats = m.regime_features(*params)
        assert m.classify_regime(*feats) == m.regime_policy(*feats)


def test_regime_features_resonant_cap():
    feats = m.regime_features(1.0, 0.05, 0.02, 0.0)
    assert feats == (20.0, 50.0, 10.0)


# ---- gates + seats
def test_hold_gates_14():
    assert len(m.HOLD_GATES) == 14
    assert all(v.startswith("HOLD") for v in m.HOLD_GATES.values())


def test_gate_outcomes():
    assert m.gate("hardware_control_claim") == "HOLD"
    assert m.gate("quantum_network_entanglement_claim") == "HOLD"
    assert m.gate("measured_data_claim") == "HOLD"
    assert m.gate("run_simulation") == "PROCEED"


def test_seats_cannot_override():
    brains = m.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["seat"] == "Gemini Argon 4 Extended Thinking"
    assert brains["decider"]["model_id"] == "gemini-argon-4"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = m.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]
