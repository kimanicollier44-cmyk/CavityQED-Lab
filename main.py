#!/usr/bin/env python3
"""CavityQED-Lab - educational cavity quantum electrodynamics lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * CavityQED-Lab simulates the Jaynes-Cummings model - one two-level
    atom coupled to one quantised cavity mode - in a TRUNCATED Fock
    basis, with every parameter a LABELLED teaching value in
    dimensionless units (the coupling g sets the scale, g = 1).
  * The dynamics are computed for real: the Hamiltonian is built from
    scratch, diagonalised by a from-scratch Jacobi eigensolver, and
    the state is evolved spectrally, so the vacuum Rabi oscillation,
    the collapse and revival, and the dispersive shifts are measured
    from the simulation and checked against the analytic formulas.
  * No device was measured, calibrated, or controlled. A cooperativity
    printed here is arithmetic on labelled rates, not a specification
    of any real cavity, atom, or circuit.

Hardware, sensor, quantum-network, single-photon-source, repeater,
fault-tolerance, and measured-data claims are permanent HOLD gates -
with Gemini Argon 4 Extended Thinking as DECIDER (names-only,
tighten-only), Muse as advisor (never overrides HOLD), and Google
Flow plan-only seats.
"""
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Labelled teaching values, dimensionless. g sets the frequency scale:
# every rate and detuning below is quoted in units of g unless a row of
# the regime table says otherwise. These are NOT measured device data.
G_DEMO = 1.0            # atom-cavity coupling at the demo point
OMEGA_C = 5.0           # labelled cavity frequency (lab-frame reports)
NBAR_DEMO = 9.0         # labelled mean photon number, coherent field
N_FOCK_VACUUM = 12      # Fock truncation for the vacuum / dressed runs
N_FOCK_COHERENT = 30    # Fock truncation for the coherent-state run
KAPPA_DEMO = 0.15       # labelled cavity decay rate for the damped run

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "real_device_calibration_claim": "HOLD - no calibrating a real cavity, atom, or circuit from labelled teaching values",
    "hardware_control_claim": "HOLD - no control of real cavity QED hardware, lasers, or cryogenic plant",
    "quantum_network_entanglement_claim": "HOLD - no claiming real remote entanglement distribution from a toy basis",
    "quantum_communication_security_claim": "HOLD - no QKD or communication-security claim from this simulation",
    "sensor_product_claim": "HOLD - no sensor product claim; a simulated cooperativity is not a sensitivity spec",
    "fault_tolerance_claim": "HOLD - no claim this toy model establishes fault-tolerant operation",
    "quantum_advantage_claim": "HOLD - no quantum-advantage claim; the whole lab runs classically in stdlib",
    "medical_imaging_claim": "HOLD - no medical imaging, sensing, or diagnostic claim",
    "single_photon_source_product_claim": "HOLD - no single-photon-source product claim from a truncated basis",
    "quantum_repeater_product_claim": "HOLD - no quantum-repeater or network-node product claim",
    "cryogenic_operation_claim": "HOLD - no operating or specifying real cryogenic hardware",
    "measured_data_claim": "HOLD - simulated values must never be presented as measured device data",
    "investment_claim": "HOLD - no investment claim; a lab model is not a prospectus",
    "publication_grade_spectroscopy_claim": "HOLD - toy spectra must not be presented as publication-grade experimental results",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Linear algebra
# The Jaynes-Cummings Hamiltonian in the |atom, n> basis is real and
# symmetric (all coupling matrix elements g*sqrt(n+1) are real), so a
# real cyclic Jacobi eigensolver is enough - and it is written from
# scratch here so every eigenvalue can be checked by reconstruction.

def jacobi_eigh(matrix, tol=1e-12, max_sweeps=100):
    """Eigen-decomposition of a real symmetric matrix, from scratch.

    Cyclic Jacobi rotations (Numerical Recipes formulation). Returns
    (eigenvalues, eigenvectors): eigenvalues sorted ascending, and
    eigenvectors[k] the matching unit vector as a plain list.
    """
    n = len(matrix)
    a = [list(row) for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = math.sqrt(sum(a[p][q] ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = a[p][p], a[q][q]
                theta = (aqq - app) / (2.0 * apq)
                t = (1.0 if theta >= 0 else -1.0) / (
                    abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                tau = s / (1.0 + c)
                a[p][p] = app - t * apq
                a[q][q] = aqq + t * apq
                a[p][q] = a[q][p] = 0.0
                for r in range(n):
                    if r == p or r == q:
                        continue
                    arp, arq = a[r][p], a[r][q]
                    a[r][p] = a[p][r] = arp - s * (arq + tau * arp)
                    a[r][q] = a[q][r] = arq + s * (arp - tau * arq)
                for r in range(n):
                    vrp, vrq = v[r][p], v[r][q]
                    v[r][p] = vrp - s * (vrq + tau * vrp)
                    v[r][q] = vrq + s * (vrp - tau * vrq)
    order = sorted(range(n), key=lambda i: a[i][i])
    eigenvalues = [a[i][i] for i in order]
    eigenvectors = [[v[r][i] for r in range(n)] for i in order]
    return eigenvalues, eigenvectors


def _golden_min(f, a, b, tol=1e-11):
    """Golden-section minimisation of a unimodal f on [a, b]."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    while abs(b - a) > tol:
        if f(c) < f(d):
            b = d
        else:
            a = c
        c = b - gr * (b - a)
        d = a + gr * (b - a)
    x = 0.5 * (a + b)
    return x, f(x)


# ================================================================ 2. Jaynes-Cummings model
# Basis: index = atom * n_fock + n, atom 0 = ground |g>, 1 = excited
# |e>, n = photon number 0..n_fock-1. The Hamiltonian conserves the
# excitation number N = a^+a + |e><e|, so it is block diagonal in
# manifolds {|e,n>, |g,n+1>} - the truncation never mixes manifolds,
# it only drops the ones above the cut.
#
# Rotating frame (omega_c = 0): H = (Delta/2) sigma_z + g (s+ a + s- a^+)
# Lab frame (omega_c > 0):      H = omega_c a^+a + (omega_a/2) sigma_z
#                                   + g (s+ a + s- a^+), omega_a = omega_c + Delta
# Populations are identical in the two frames; the lab frame is used
# only where an absolute dressed energy is the point of the report.

def build_jc_hamiltonian(g, delta, n_fock, omega_c=0.0):
    """Full Jaynes-Cummings Hamiltonian in the truncated Fock basis."""
    dim = 2 * n_fock
    h = [[0.0] * dim for _ in range(dim)]
    omega_a = omega_c + delta
    for n in range(n_fock):
        h[n][n] = omega_c * n - 0.5 * omega_a                    # |g, n>
        h[n_fock + n][n_fock + n] = omega_c * n + 0.5 * omega_a  # |e, n>
    for n in range(n_fock - 1):
        e_n = n_fock + n          # |e, n>
        g_n1 = n + 1              # |g, n+1>
        h[e_n][g_n1] = h[g_n1][e_n] = g * math.sqrt(n + 1)
    return h


def manifold_eigenvalues(g, delta, n):
    """Exact eigenvalues of one manifold {|e,n>, |g,n+1>} (rotating frame)."""
    half = 0.5 * math.sqrt(delta * delta + 4.0 * g * g * (n + 1))
    return (-half, half)


def _pe_function(eigenvalues, eigenvectors, psi0, n_fock):
    """Excited-state probability P_e(t) as an exact cosine sum.

    P_e(t) = sum_{k,l} a_k a_l S^e_{kl} cos((lam_k - lam_l) t) with
    a_k = <v_k|psi0> and S^e_{kl} the excited-atom overlap of the two
    eigenvectors - so each evaluation is real arithmetic, no complex
    exponentials, and the result is deterministic by construction.
    """
    dim = len(psi0)
    amps = [sum(eigenvectors[k][j] * psi0[j] for j in range(dim))
            for k in range(dim)]
    terms = []
    for k in range(dim):
        for l in range(k, dim):
            s_e = sum(eigenvectors[k][n_fock + n] * eigenvectors[l][n_fock + n]
                      for n in range(n_fock))
            coef = amps[k] * amps[l] * s_e
            if abs(coef) < 1e-15:
                continue
            terms.append((eigenvalues[k] - eigenvalues[l], coef,
                          1.0 if k == l else 2.0))

    def pe(t):
        return sum(mult * coef * math.cos(freq * t)
                   for freq, coef, mult in terms)

    return pe


def vacuum_state(n_fock):
    """|e, 0>: the atom excited, the cavity empty."""
    psi = [0.0] * (2 * n_fock)
    psi[n_fock] = 1.0
    return psi


def coherent_state(nbar, n_fock):
    """|g> x |alpha>, alpha = sqrt(nbar) real: Poisson photon weights."""
    psi = [0.0] * (2 * n_fock)
    for n in range(n_fock):
        psi[n] = math.exp(-0.5 * nbar) * nbar ** (0.5 * n) / math.sqrt(
            math.factorial(n))
    return psi


# ================================================================ 3. Vacuum Rabi + dressed states
def vacuum_rabi_report(g=G_DEMO, n_fock=N_FOCK_VACUUM):
    """Vacuum Rabi oscillation measured from the simulated P_e(t).

    The frequency is MEASURED: the first two minima of P_e(t) are
    located by golden-section search and the period is their distance.
    For the resonant JC model P_e(t) = cos^2(g t), so the measured
    angular frequency must land on 2g - the checks assert that it does.
    """
    h = build_jc_hamiltonian(g, 0.0, n_fock)
    eigenvalues, eigenvectors = jacobi_eigh(h)
    pe = _pe_function(eigenvalues, eigenvectors, vacuum_state(n_fock),
                      n_fock)
    t_min1, p_min1 = _golden_min(pe, 0.5 / g, 2.5 / g)
    t_min2, _ = _golden_min(pe, 3.5 / g, 5.5 / g)
    period = t_min2 - t_min1
    return {
        "g": g,
        "n_fock": n_fock,
        "p_e_at_0": round(pe(0.0), 12),
        "first_min_time": round(t_min1, 6),
        "first_min_value": round(p_min1, 12),
        "p_e_at_half_period": round(pe(math.pi / g), 12),
        "period_measured": round(period, 6),
        "frequency_measured": round(2.0 * math.pi / period, 6),
        "frequency_analytic_2g": round(2.0 * g, 6),
    }


def dressed_splitting_report(g=G_DEMO, delta=0.0, omega_c=OMEGA_C,
                             n_fock=N_FOCK_VACUUM):
    """Vacuum Rabi splitting from the lab-frame dressed eigenvalues.

    The N = 1 doublet {|e,0>, |g,1>} sits at centre omega_c/2; the two
    eigenvalues nearest the centre are the dressed pair and their
    distance is the splitting - sqrt(Delta^2 + 4 g^2), i.e. 2g on
    resonance. Nothing is quoted: the eigenvalues come out of the
    Jacobi solver on the full truncated Hamiltonian.
    """
    h = build_jc_hamiltonian(g, delta, n_fock, omega_c=omega_c)
    eigenvalues, _ = jacobi_eigh(h)
    centre = 0.5 * omega_c
    pair = sorted(eigenvalues, key=lambda e: abs(e - centre))[:2]
    splitting = abs(pair[1] - pair[0])
    return {
        "g": g,
        "delta": delta,
        "omega_c": omega_c,
        "doublet_centre": centre,
        "dressed_pair": [round(e, 6) for e in sorted(pair)],
        "splitting_measured": round(splitting, 6),
        "splitting_analytic": round(math.sqrt(delta * delta + 4.0 * g * g), 6),
    }


# ================================================================ 4. Collapse and revival
_REVIVAL_CACHE = {}


def revival_report(nbar=NBAR_DEMO, g=G_DEMO, n_fock=N_FOCK_COHERENT):
    """Collapse and revival with a coherent-state field.

    The atom starts in |g>, the field in |alpha> with labelled mean
    photon number nbar (Poisson weights). Each photon-number component
    Rabi-oscillates at its own frequency 2g*sqrt(n+1); the components
    dephase (collapse) and rephase (revival). The revival time is
    MEASURED as the highest P_e peak inside [0.7, 1.3] x t_estimate
    and compared with the analytic estimate 2*pi*sqrt(nbar)/g.
    """
    key = (nbar, g, n_fock)
    if key in _REVIVAL_CACHE:
        return _REVIVAL_CACHE[key]
    h = build_jc_hamiltonian(g, 0.0, n_fock)
    eigenvalues, eigenvectors = jacobi_eigh(h)
    psi0 = coherent_state(nbar, n_fock)
    pe = _pe_function(eigenvalues, eigenvectors, psi0, n_fock)
    t_est = 2.0 * math.pi * math.sqrt(nbar) / g
    # Early oscillation amplitude (first lobe).
    early = [pe(0.05 * i) for i in range(41)]            # t in [0, 2]
    # Collapse window: statistics of P_e where the components dephase.
    lo, hi = 0.25 * t_est, 0.55 * t_est
    window = [pe(lo + (hi - lo) * i / 200.0) for i in range(201)]
    # Revival: grid argmax, then golden-section refinement.
    grid = [0.7 * t_est + 0.6 * t_est * i / 240.0 for i in range(241)]
    t_grid = max(grid, key=pe)
    t_rev, _ = _golden_min(lambda t: -pe(t), t_grid - 0.3, t_grid + 0.3)
    p_rev = pe(t_rev)
    report = {
        "nbar_labelled": nbar,
        "g": g,
        "n_fock": n_fock,
        "poisson_norm_truncated": round(sum(x * x for x in psi0), 12),
        "early_max": round(max(early), 6),
        "early_min": round(min(early), 6),
        "collapse_mean": round(sum(window) / len(window), 6),
        "collapse_min": round(min(window), 6),
        "collapse_max": round(max(window), 6),
        "collapse_range": round(max(window) - min(window), 6),
        "revival_time_estimate": round(t_est, 6),
        "revival_time_measured": round(t_rev, 6),
        "revival_peak": round(p_rev, 6),
        "revival_time_rel_error": round(abs(t_rev - t_est) / t_est, 6),
    }
    _REVIVAL_CACHE[key] = report
    return report


# ================================================================ 5. Coupling regimes
# Labelled parameter table (dimensionless, g-scale). Every rate here
# is a teaching value chosen to put the row in its named regime - the
# table demonstrates the arithmetic of the regime boundaries, it does
# not describe any real device.
REGIME_TABLE = [
    {"name": "demo_point", "g": 1.0, "kappa": 0.05, "gamma": 0.02,
     "delta": 0.0},
    {"name": "high_cooperativity_labelled", "g": 0.8, "kappa": 0.02,
     "gamma": 0.01, "delta": 0.0},
    {"name": "marginal_labelled", "g": 0.3, "kappa": 0.3, "gamma": 0.1,
     "delta": 0.0},
    {"name": "lossy_cavity_labelled", "g": 0.2, "kappa": 0.5,
     "gamma": 0.1, "delta": 0.0},
    {"name": "bad_cavity_labelled", "g": 0.1, "kappa": 2.0,
     "gamma": 0.05, "delta": 0.0},
]


def cooperativity(g, kappa, gamma):
    """C = g^2 / (kappa * gamma): coherent coupling vs both loss rates."""
    return g * g / (kappa * gamma)


def cavity_induced_decay(g, kappa):
    """Bad-cavity-limit cavity-induced decay rate 4 g^2 / kappa.

    Valid only when kappa >> g (the cavity empties before one Rabi
    cycle completes) - the label is part of the definition here.
    """
    return 4.0 * g * g / kappa


def purcell_factor(g, kappa, gamma):
    """Rate-model Purcell factor F_p = Gamma_cav / gamma = 4C.

    The enhancement of the atomic emission rate into the cavity mode
    over the free-space rate gamma, in the bad-cavity rate model.
    """
    return cavity_induced_decay(g, kappa) / gamma


def strong_coupling(g, kappa, gamma):
    """Strong coupling: one coherent exchange beats both loss channels."""
    return g > kappa and g > gamma


def splitting_resolved(g, kappa, gamma):
    """The vacuum Rabi doublet is resolvable: 2g > (kappa + gamma)/2."""
    return 2.0 * g > 0.5 * (kappa + gamma)


def regimes_report():
    rows = []
    for row in REGIME_TABLE:
        g, kappa, gamma = row["g"], row["kappa"], row["gamma"]
        rows.append({
            **row,
            "cooperativity": round(cooperativity(g, kappa, gamma), 6),
            "purcell_factor": round(purcell_factor(g, kappa, gamma), 6),
            "cavity_induced_decay": round(cavity_induced_decay(g, kappa), 6),
            "strong_coupling": strong_coupling(g, kappa, gamma),
            "splitting_resolved": splitting_resolved(g, kappa, gamma),
        })
    return {
        "rows": rows,
        "cooperativity_max": max(r["cooperativity"] for r in rows),
        "cooperativity_min": min(r["cooperativity"] for r in rows),
        "note": "labelled dimensionless teaching values - not device data",
    }


# ================================================================ 6. Dispersive regime
def dispersive_report(g=G_DEMO, delta=10.0, n_fock=8):
    """Dispersive regime: |Delta| >> g, no energy exchange, only pulls.

    chi = g^2 / Delta is the cavity pull per atomic state; the atomic
    line is number-split by 2*chi per photon (the ac Stark shift per
    photon in this convention). The analytic shifts are checked
    against the EXACT manifold eigenvalues of the same model: the
    atom-like dressed shift for manifold n is
    (1/2)(sqrt(Delta^2 + 4 g^2 (n+1)) - Delta), and its difference
    between neighbouring manifolds is the measured per-photon pull.
    The full truncated spectrum is also compared, eigenvalue by
    eigenvalue, against the union of all manifold pairs plus the two
    uncoupled edge states |g,0> and |e, n_fock-1>.
    """
    chi = g * g / delta
    shifts = []
    for n in range(3):
        exact = manifold_eigenvalues(g, delta, n)[1] - 0.5 * delta
        shifts.append({
            "n": n,
            "dressed_shift_exact": round(exact, 6),
            "dressed_shift_dispersive": round(chi * (n + 1), 6),
        })
    per_photon_exact = shifts[1]["dressed_shift_exact"] \
        - shifts[0]["dressed_shift_exact"]
    h = build_jc_hamiltonian(g, delta, n_fock)
    eigenvalues, _ = jacobi_eigh(h)
    expected = [-0.5 * delta, 0.5 * delta]
    for n in range(n_fock - 1):
        expected.extend(manifold_eigenvalues(g, delta, n))
    expected.sort()
    max_dev = max(abs(a - b) for a, b in zip(sorted(eigenvalues), expected))
    return {
        "g": g,
        "delta": delta,
        "validity_g_over_delta": round(g / delta, 6),
        "chi": round(chi, 6),
        "cavity_pull_per_state": round(chi, 6),
        "number_splitting_2chi": round(2.0 * chi, 6),
        "ac_stark_per_photon": round(2.0 * chi, 6),
        "manifold_shifts": shifts,
        "per_photon_pull_exact": round(per_photon_exact, 6),
        "full_spectrum_max_deviation": max_dev,
    }


# ================================================================ 7. Dissipation (Lindblad-lite)
# A deliberately small open-system model: the three states the vacuum
# dynamics can reach, {|e,0>, |g,1>, |g,0>}, a density matrix, and two
# collapse channels (cavity decay kappa: |g,1> -> |g,0>; atomic decay
# gamma: |e,0> -> |g,0>), integrated by fixed-step RK4. "Lindblad-lite"
# means exactly this: one Markovian channel per loss, no bath memory,
# no drive - enough to watch kappa damp the Rabi oscillation, labelled
# as the simplification it is.

def _cmat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]


def _cmat_dagger(a):
    return [[a[j][i].conjugate() for j in range(3)] for i in range(3)]


def _lindblad_rhs(rho, h, channels):
    comm = _cmat_mul(h, rho)
    comm2 = _cmat_mul(rho, h)
    out = [[-1j * (comm[i][j] - comm2[i][j]) for j in range(3)]
           for i in range(3)]
    for rate, c in channels:
        if rate == 0.0:
            continue
        cd = _cmat_dagger(c)
        cdc = _cmat_mul(cd, c)
        jump = _cmat_mul(_cmat_mul(c, rho), cd)
        anti1 = _cmat_mul(cdc, rho)
        anti2 = _cmat_mul(rho, cdc)
        for i in range(3):
            for j in range(3):
                out[i][j] += rate * (jump[i][j]
                                     - 0.5 * (anti1[i][j] + anti2[i][j]))
    return out


def damped_rabi_report(g=G_DEMO, kappa=KAPPA_DEMO, gamma=0.0,
                       t_max=40.0, step=0.01):
    """Damped vacuum Rabi oscillation under cavity decay kappa.

    Returns the successive P_e peak heights (the measured envelope),
    a log-linear fit of the envelope decay rate - expected near
    kappa/2, since the excitation spends about half its time in the
    lossy cavity - and the worst trace drift of the integration.
    """
    zero = [[0j] * 3 for _ in range(3)]
    h = [[0j] * 3 for _ in range(3)]
    h[0][1] = h[1][0] = g
    cav = [row[:] for row in zero]
    cav[2][1] = 1.0            # a: |g,1> -> |g,0>
    atom = [row[:] for row in zero]
    atom[2][0] = 1.0           # sigma-: |e,0> -> |g,0>
    channels = [(kappa, cav), (gamma, atom)]
    rho = [[0j] * 3 for _ in range(3)]
    rho[0][0] = 1.0 + 0j       # start in |e,0>
    steps = int(round(t_max / step))
    series = [(0.0, 1.0)]
    max_trace_err = 0.0
    for s in range(1, steps + 1):
        k1 = _lindblad_rhs(rho, h, channels)
        r2 = [[rho[i][j] + 0.5 * step * k1[i][j] for j in range(3)]
              for i in range(3)]
        k2 = _lindblad_rhs(r2, h, channels)
        r3 = [[rho[i][j] + 0.5 * step * k2[i][j] for j in range(3)]
              for i in range(3)]
        k3 = _lindblad_rhs(r3, h, channels)
        r4 = [[rho[i][j] + step * k3[i][j] for j in range(3)]
              for i in range(3)]
        k4 = _lindblad_rhs(r4, h, channels)
        rho = [[rho[i][j] + step / 6.0 * (k1[i][j] + 2 * k2[i][j]
                                          + 2 * k3[i][j] + k4[i][j])
                for j in range(3)] for i in range(3)]
        trace = sum(rho[i][i] for i in range(3)).real
        max_trace_err = max(max_trace_err, abs(trace - 1.0))
        series.append((s * step, rho[0][0].real))
    peaks = []
    for i in range(1, len(series) - 1):
        if series[i][1] >= series[i - 1][1] and series[i][1] > series[i + 1][1] \
                and series[i][0] > 0.5:
            peaks.append((round(series[i][0], 4), round(series[i][1], 6)))
    fit_rate = None
    if len(peaks) >= 3:
        xs = [p[0] for p in peaks]
        ys = [math.log(p[1]) for p in peaks]
        xm = sum(xs) / len(xs)
        ym = sum(ys) / len(ys)
        slope = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) \
            / sum((x - xm) ** 2 for x in xs)
        fit_rate = round(-slope, 6)
    return {
        "g": g,
        "kappa_labelled": kappa,
        "gamma_labelled": gamma,
        "peaks": peaks,
        "peak_count": len(peaks),
        "envelope_decay_rate_fit": fit_rate,
        "envelope_decay_rate_expected_kappa_over_2": round(kappa / 2.0, 6),
        "final_p_e": round(series[-1][1], 6),
        "max_trace_error": max_trace_err,
    }


# ================================================================ 8. kNN regime classifier
# Features: (g/kappa, g/gamma, g/|Delta|) - the three ratios that
# decide the regime. The training grid is labelled by the lab's own
# regime policy (regime_policy below): the classifier interpolates
# the policy, it does not replace the boundary arithmetic of section 5.

def regime_policy(g_over_kappa, g_over_gamma, g_over_delta):
    """The lab's regime policy, distilled to four labels.

    Order matters and is stated openly: a cavity that empties faster
    than one Rabi cycle (g/kappa < 1) on a near-resonant atom is the
    bad-cavity (Purcell) regime; a far-detuned atom (g/|Delta| <= 0.2)
    is dispersive; coupling beating both losses is strong coupling;
    anything else is weak coupling.
    """
    if g_over_kappa < 1.0 and g_over_delta > 0.2:
        return "bad-cavity"
    if g_over_delta <= 0.2:
        return "dispersive"
    if g_over_kappa >= 1.0 and g_over_gamma >= 1.0:
        return "strong-coupling"
    return "weak-coupling"


def regime_features(g, kappa, gamma, delta):
    """Regime ratios for one parameter set (Delta = 0 -> ratio capped)."""
    z = g / abs(delta) if delta != 0 else 10.0
    return (g / kappa, g / gamma, min(z, 10.0))


def _training_grid():
    rows = []
    for x in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 25.0):
        for y in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 25.0):
            for z in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 10.0):
                rows.append(((x, y, z), regime_policy(x, y, z)))
    return rows


TRAINING_GRID = _training_grid()


def _scaled(features):
    x, y, z = features
    return (math.log10(x) / 2.0, math.log10(y) / 2.0, math.log10(z) / 2.0)


def classify_regime(g_over_kappa, g_over_gamma, g_over_delta, k=5):
    """kNN over the policy-labelled grid, deterministic tie-breaking."""
    target = _scaled((g_over_kappa, g_over_gamma, g_over_delta))
    ranked = []
    for features, label in TRAINING_GRID:
        point = _scaled(features)
        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(target, point)))
        ranked.append((dist, features, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for _, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0) + 1
    priority = {"strong-coupling": 0, "dispersive": 1,
                "weak-coupling": 2, "bad-cavity": 3}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


# The lab's own named parameter sets, classified.
CLASSIFIER_SCENARIOS = {
    "demo_point": (1.0, 0.05, 0.02, 0.0),
    "dispersive_point": (1.0, 0.05, 0.02, 10.0),
    "bad_cavity_point": (0.1, 2.0, 0.05, 0.0),
    "weak_point": (0.1, 0.05, 0.5, 0.1),
}


def classifier_report():
    report = {}
    for name, (g, kappa, gamma, delta) in CLASSIFIER_SCENARIOS.items():
        features = regime_features(g, kappa, gamma, delta)
        report[name] = {
            "params": {"g": g, "kappa": kappa, "gamma": gamma,
                       "delta": delta},
            "features": [round(f, 6) for f in features],
            "label": classify_regime(*features),
        }
    return report


# ================================================================ 9. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="cavityqed-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "role": "DECIDER",
        "model_id": "gemini-argon-4",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all hardware and device gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="cavityqed-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Say it in the first line of every report: every rate is a labelled dimensionless teaching value and the basis is truncated - then show the measured Rabi frequency next to 2g and the revival time next to its estimate, so the simulation grades itself.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["two cavity mirrors facing each other across a dark gap",
                  "one atom between the mirrors trading a single photon back and forth",
                  "a wave of many photon numbers smearing out, then rephasing into one pulse",
                  "the cavity line splitting in two as the atom crosses resonance",
                  "a damped oscillation spiral settling as the cavity leaks"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 1. eigensolver + Hamiltonian construction
    evals, _ = jacobi_eigh([[2.0, 0.0], [0.0, -3.0]])
    chk("jacobi_diagonal", abs(evals[0] + 3.0) < 1e-12
        and abs(evals[1] - 2.0) < 1e-12)
    h = build_jc_hamiltonian(G_DEMO, 0.0, N_FOCK_VACUUM)
    evals, evecs = jacobi_eigh(h)
    dim = 2 * N_FOCK_VACUUM
    recon = max(abs(sum(evecs[k][i] * evals[k] * evecs[k][j]
                        for k in range(dim)) - h[i][j])
                for i in range(dim) for j in range(dim))
    chk("jacobi_reconstruction", recon < 1e-9)
    ortho = max(abs(sum(evecs[k][i] * evecs[l][i] for i in range(dim))
                    - (1.0 if k == l else 0.0))
                for k in range(dim) for l in range(dim))
    chk("jacobi_orthonormal", ortho < 1e-9)
    chk("jc_coupling_element",
        h[N_FOCK_VACUUM][1] == G_DEMO * math.sqrt(1)
        and h[N_FOCK_VACUUM + 1][2] == G_DEMO * math.sqrt(2))
    chk("jc_excitation_conserved",
        h[N_FOCK_VACUUM][2] == 0.0 and h[0][N_FOCK_VACUUM + 1] == 0.0
        and h[3][3] == 0.0)

    # 2. vacuum Rabi - measured against 2g
    rabi = vacuum_rabi_report()
    chk("rabi_starts_excited", rabi["p_e_at_0"] == 1.0)
    chk("rabi_first_min_is_zero", rabi["first_min_value"] < 1e-9)
    chk("rabi_first_min_time", abs(rabi["first_min_time"]
                                   - math.pi / (2 * G_DEMO)) < 1e-4)
    chk("rabi_full_period_return", rabi["p_e_at_half_period"] > 1 - 1e-9)
    chk("rabi_frequency_is_2g",
        abs(rabi["frequency_measured"] - 2 * G_DEMO) / (2 * G_DEMO) < 0.005)
    chk("rabi_period_is_pi_over_g",
        abs(rabi["period_measured"] - math.pi / G_DEMO) < 1e-3)

    # 3. dressed splitting from lab-frame eigenvalues
    dress0 = dressed_splitting_report(delta=0.0)
    chk("splitting_resonant_is_2g",
        abs(dress0["splitting_measured"] - 2 * G_DEMO) < 1e-6)
    dress3 = dressed_splitting_report(delta=3.0)
    chk("splitting_detuned_formula",
        abs(dress3["splitting_measured"] - math.sqrt(13.0)) < 1e-6
        and abs(dress3["splitting_measured"]
                - dress3["splitting_analytic"]) < 1e-9)

    # 4. collapse and revival
    rev = revival_report()
    chk("poisson_norm_kept", rev["poisson_norm_truncated"] > 1 - 1e-7)
    chk("revival_early_oscillation",
        rev["early_max"] > 0.9 and rev["early_min"] < 0.1)
    chk("collapse_range_small", rev["collapse_range"] < 0.15)
    chk("collapse_below_early_amplitude",
        rev["collapse_max"] < rev["early_max"] - 0.2)
    chk("revival_peak_high", rev["revival_peak"] > 0.75)
    chk("revival_time_matches_estimate", rev["revival_time_rel_error"] < 0.12)

    # 5. regimes: cooperativity, Purcell, boundaries
    reg = regimes_report()
    rows = {r["name"]: r for r in reg["rows"]}
    chk("cooperativity_demo", rows["demo_point"]["cooperativity"] == 1000.0)
    chk("cooperativity_extremes",
        reg["cooperativity_max"] == 3200.0
        and reg["cooperativity_min"] == 0.1)
    chk("purcell_is_4C",
        all(abs(r["purcell_factor"] - 4 * r["cooperativity"]) < 1e-9
            for r in reg["rows"]))
    chk("purcell_bad_cavity",
        rows["bad_cavity_labelled"]["purcell_factor"] == 0.4
        and rows["bad_cavity_labelled"]["cavity_induced_decay"] == 0.02)
    chk("strong_coupling_flags",
        rows["demo_point"]["strong_coupling"] is True
        and rows["bad_cavity_labelled"]["strong_coupling"] is False
        and rows["lossy_cavity_labelled"]["strong_coupling"] is False)
    chk("splitting_resolved_flags",
        rows["demo_point"]["splitting_resolved"] is True
        and rows["marginal_labelled"]["splitting_resolved"] is True)

    # 6. dispersive regime
    disp = dispersive_report()
    chk("chi_value", disp["chi"] == 0.1
        and disp["number_splitting_2chi"] == 0.2
        and disp["ac_stark_per_photon"] == 0.2)
    chk("dispersive_validity", disp["validity_g_over_delta"] == 0.1
        and disp["validity_g_over_delta"] < 0.2)
    s0 = disp["manifold_shifts"][0]
    chk("dispersive_shift_n0",
        abs(s0["dressed_shift_exact"] - s0["dressed_shift_dispersive"])
        / s0["dressed_shift_dispersive"] < 0.05)
    chk("dispersive_per_photon_pull",
        abs(disp["per_photon_pull_exact"] - disp["chi"]) / disp["chi"] < 0.05)
    chk("full_spectrum_matches_manifolds",
        disp["full_spectrum_max_deviation"] < 1e-8)

    # 7. dissipation
    damped = damped_rabi_report()
    chk("lindblad_trace_conserved", damped["max_trace_error"] < 1e-9)
    heights = [p[1] for p in damped["peaks"]]
    chk("damped_peaks_decrease",
        len(heights) >= 5
        and all(heights[i + 1] < heights[i] for i in range(len(heights) - 1)))
    chk("damped_envelope_rate",
        damped["envelope_decay_rate_fit"] is not None
        and 0.5 * KAPPA_DEMO / 2 < damped["envelope_decay_rate_fit"]
        < 1.6 * KAPPA_DEMO / 2)
    chk("damped_final_small", damped["final_p_e"] < 0.05)
    undamped = damped_rabi_report(kappa=0.0, t_max=10.0)
    chk("undamped_peaks_stay_one",
        all(abs(p[1] - 1.0) < 1e-4 for p in undamped["peaks"])
        and undamped["peak_count"] >= 3)

    # 8. classifier - labels come from the regime policy itself
    rep = classifier_report()
    chk("classifier_demo_strong",
        rep["demo_point"]["label"] == "strong-coupling")
    chk("classifier_dispersive",
        rep["dispersive_point"]["label"] == "dispersive")
    chk("classifier_bad_cavity",
        rep["bad_cavity_point"]["label"] == "bad-cavity")
    chk("classifier_weak", rep["weak_point"]["label"] == "weak-coupling")
    chk("classifier_matches_policy",
        all(classify_regime(*regime_features(*params))
            == regime_policy(*regime_features(*params))
            for params in CLASSIFIER_SCENARIOS.values()))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_hardware_hold", gate("hardware_control_claim") == "HOLD")
    chk("gate_network_hold",
        gate("quantum_network_entanglement_claim") == "HOLD")
    chk("gate_measured_data_hold", gate("measured_data_claim") == "HOLD")
    chk("gate_sensor_hold", gate("sensor_product_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking"
        and brains["decider"]["model_id"] == "gemini-argon-4")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("CavityQED-Lab - cavity quantum electrodynamics lab")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - a truncated toy, graded against analytics")
    print("    One atom, one cavity mode, the Jaynes-Cummings Hamiltonian in")
    print("    a truncated Fock basis. Every rate is a LABELLED dimensionless")
    print("    teaching value (g = 1 sets the scale). No device was measured")
    print("    or calibrated; the simulation is graded against the analytic")
    print("    formulas it is supposed to reproduce, in the open.")
    rabi = vacuum_rabi_report()
    dress = dressed_splitting_report()
    print("\n[2] Vacuum Rabi oscillation - measured from the simulation")
    print(f"    P_e(0) = {rabi['p_e_at_0']:.6f} | first minimum "
          f"{rabi['first_min_value']:.2e} at t = {rabi['first_min_time']:.4f}")
    print(f"    period measured {rabi['period_measured']:.4f} "
          f"(pi/g = {math.pi / G_DEMO:.4f})")
    print(f"    frequency measured {rabi['frequency_measured']:.4f} "
          f"vs analytic 2g = {rabi['frequency_analytic_2g']:.4f}")
    print(f"    dressed doublet (lab frame, omega_c = {OMEGA_C}): "
          f"{dress['dressed_pair']} -> splitting "
          f"{dress['splitting_measured']:.4f} = 2g")
    rev = revival_report()
    print("\n[3] Collapse and revival - coherent field, labelled nbar = "
          f"{NBAR_DEMO:.0f}")
    print(f"    early oscillation: max {rev['early_max']:.4f} / "
          f"min {rev['early_min']:.4f}")
    print(f"    collapse window: mean {rev['collapse_mean']:.4f}, "
          f"range {rev['collapse_range']:.4f} "
          f"[{rev['collapse_min']:.4f}, {rev['collapse_max']:.4f}]")
    print(f"    revival peak {rev['revival_peak']:.4f} at "
          f"t = {rev['revival_time_measured']:.4f} | estimate "
          f"2*pi*sqrt(nbar)/g = {rev['revival_time_estimate']:.4f} "
          f"(rel err {rev['revival_time_rel_error']:.4f})")
    reg = regimes_report()
    print("\n[4] Coupling regimes - cooperativity on labelled rates")
    for row in reg["rows"]:
        print(f"    {row['name']:28s} C = {row['cooperativity']:9.2f} | "
              f"Purcell {row['purcell_factor']:9.2f} | strong "
              f"{row['strong_coupling']} | resolved "
              f"{row['splitting_resolved']}")
    print(f"    table extremes: C from {reg['cooperativity_min']} to "
          f"{reg['cooperativity_max']} (arithmetic on labels, not devices)")
    disp = dispersive_report()
    print("\n[5] Dispersive regime - g/Delta = "
          f"{disp['validity_g_over_delta']:.2f} << 1")
    print(f"    chi = g^2/Delta = {disp['chi']:.4f} | number splitting "
          f"2*chi = {disp['number_splitting_2chi']:.4f} | ac Stark per "
          f"photon {disp['ac_stark_per_photon']:.4f}")
    print(f"    exact per-photon pull {disp['per_photon_pull_exact']:.4f} "
          f"vs chi {disp['chi']:.4f} | full spectrum vs manifolds: "
          f"max dev {disp['full_spectrum_max_deviation']:.2e}")
    damped = damped_rabi_report()
    print(f"\n[6] Dissipation - Lindblad-lite, labelled kappa = {KAPPA_DEMO}")
    print(f"    P_e peaks: {[(t, round(p, 4)) for t, p in damped['peaks'][:5]]}")
    print(f"    envelope decay fit {damped['envelope_decay_rate_fit']} "
          f"vs kappa/2 = {damped['envelope_decay_rate_expected_kappa_over_2']} "
          f"| final P_e {damped['final_p_e']} | trace err "
          f"{damped['max_trace_error']:.2e}")
    print("\n[7] Regime classifier (TinyML kNN, labels from the lab policy)")
    for name, entry in classifier_report().items():
        print(f"    {name:36s} -> {entry['label']}")
    print("\n[8] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[9] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "rabi":
        print(json.dumps({"vacuum_rabi": vacuum_rabi_report(),
                          "dressed_splitting": dressed_splitting_report()},
                         indent=2, sort_keys=True))
        return 0
    if cmd == "revival":
        print(json.dumps(revival_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "regimes":
        print(json.dumps(regimes_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "dispersive":
        print(json.dumps(dispersive_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "decay":
        print(json.dumps(damped_rabi_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "hardware_control_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
