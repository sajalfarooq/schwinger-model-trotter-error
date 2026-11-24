"""Compare ideal product formulas under a common two-qubit rotation budget."""
import numpy as np
import schwinger as pilot


def schedule(gates, steps, symmetric):
    half = [(s, p, c / 2) for s, p, c in gates]
    sequence = (half + half[::-1] if symmetric else gates) * steps
    merged = []
    for sites, pauli, coefficient in sequence:
        if merged and merged[-1][1] is pauli:
            merged[-1][2] += coefficient
        else:
            merged.append([sites, pauli, coefficient])
    return merged


def run(n):
    pilot.N = n
    h, gates, vacuum, obs = pilot.model()
    times = np.linspace(0, 2, 81)
    energy, basis = np.linalg.eigh(h)
    exact = (np.exp(-1j * times[:, None] * energy) * (basis.conj().T @ vacuum)) @ basis.T
    reference = abs(exact)**2 @ obs.T
    rows = []
    # Each Pauli string has one nonzero entry per row.
    actions = {}
    for _, p, _ in gates:
        indices = np.argmax(abs(p), axis=1)
        actions[id(p)] = indices, p[np.arange(len(p)), indices]
    for symmetric in (False, True):
        for steps in (2, 4, 8, 16, 32):
            circuit = schedule(gates, steps, symmetric)
            psi = np.tile(vacuum, (len(times), 1)).astype(complex)
            for sites, p, c in circuit:
                indices, phase = actions[id(p)]
                angle = times[:, None] * c / steps
                psi = np.cos(angle)*psi - 1j*np.sin(angle)*phase*psi[:, indices]
            assert np.allclose(np.sum(abs(psi)**2, axis=1), 1)
            error = abs(abs(psi)**2 @ obs.T - reference)
            rows.append((n, 2 if symmetric else 1, steps,
                         sum(len(s)==2 for s, _, _ in circuit),
                         error[:, 0].max(), error[:, 1:].max()))
    return rows


if __name__ == "__main__":
    rows = run(4) + run(6)
    pilot.write_csv("formula_comparison.csv",
                    ["sites","order","steps","two_qubit_rotations","max_density_error","max_flux_error"], rows)
    for row in rows:
        print("%d %d %d %d %.6f %.6f" % row)
