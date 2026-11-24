"""Four-site Schwinger pilot using NumPy and a native Pauli rotation gate set."""
import csv
from pathlib import Path

import numpy as np

N = 4
W, MASS, J, BACKGROUND = 1.0, 0.5, 1.0, 0.0
PAULI = dict(zip("IXYZ", map(np.array, [
    [[1, 0], [0, 1]], [[0, 1], [1, 0]],
    [[0, -1j], [1j, 0]], [[1, 0], [0, -1]],
])))


def operator(sites):
    result = np.ones((1, 1))
    for site in range(N):
        result = np.kron(result, PAULI[sites.get(site, "I")])
    return result


def model():
    # One-based staggering as in Martinez et al. (+Z on site 1 in vacuum).
    stagger = (-1) ** np.arange(1, N + 1)
    spins = np.array([operator({site: "Z"}).diagonal() for site in range(N)])
    charges = (spins + stagger[:, None]) / 2
    fields = BACKGROUND + charges.cumsum(axis=0)[:-1]
    density = (1 + stagger[:, None] * spins).mean(axis=0) / 2
    diagonal = MASS / 2 * (stagger @ spins) + J * (fields**2).sum(axis=0)
    hamiltonian = np.diag(diagonal)
    gates = []
    supports = [(i,) for i in range(N)]
    supports += [(i, j) for i in range(N) for j in range(i + 1, N)]
    for sites in supports:
        pauli = operator(dict.fromkeys(sites, "Z"))
        coefficient = np.mean(pauli.diagonal() * diagonal)
        if abs(coefficient) > 1e-12:
            gates.append((sites, pauli, coefficient))
    offset = diagonal.mean()
    assert np.allclose(sum(c * p for _, p, c in gates) + offset * np.eye(2**N),
                       hamiltonian)
    hamiltonian = hamiltonian.astype(complex)
    # Build benchmark hopping independently by swapping unlike neighboring bits.
    for i in range(N - 1):
        mask = (1 << (N - 1 - i)) | (1 << (N - 2 - i))
        for state in range(2**N):
            if (state & mask).bit_count() == 1:
                hamiltonian[state ^ mask, state] += W
        for axis in "XY":
            gates.append(((i, i + 1), operator({i: axis, i + 1: axis}), W / 2))
    vacuum = np.eye(2**N)[:, np.flatnonzero(density == 0).item()]
    assert np.allclose(fields[:, density == 0], BACKGROUND)
    total_charge = charges.sum(axis=0)
    assert np.allclose(hamiltonian * (total_charge[:, None] - total_charge), 0)
    return hamiltonian, gates, vacuum, np.vstack((density, fields))


def circuit(gates, vacuum, time, steps, noise):
    rho = np.outer(vacuum, vacuum).astype(complex)
    layers = [0] * N
    two_qubit_gates = 0
    identity = np.eye(2**N)
    rotations = []
    for sites, pauli, coefficient in gates:
        angle = coefficient * time / steps
        unitary = np.cos(angle) * identity - 1j * np.sin(angle) * pauli
        errors = []
        if noise and len(sites) == 2:
            errors = [operator(dict(zip(sites, (a, b)))) for a in PAULI for b in PAULI]
        rotations.append((sites, unitary, unitary.conj().T, errors))

    for _ in range(steps):
        for sites, unitary, adjoint, errors in rotations:
            rho = unitary @ rho @ adjoint
            if errors:
                # Uniform twirl, including identity, replaces the pair by I/4.
                mixed = sum(error @ rho @ error for error in errors) / 16
                rho = (1 - noise) * rho + noise * mixed
            layer = 1 + max(layers[site] for site in sites)
            for site in sites:
                layers[site] = layer
            two_qubit_gates += len(sites) == 2
    assert np.isclose(np.trace(rho), 1)
    return rho, max(layers), two_qubit_gates


def write_csv(name, columns, rows):
    with Path(__file__).with_name(name).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    hamiltonian, gates, vacuum, observables = model()
    energy, eigenvectors = np.linalg.eigh(hamiltonian)
    times = np.linspace(0, 2, 21)
    states = (np.exp(-1j * times[:, None] * energy) *
              (eigenvectors.conj().T @ vacuum)) @ eigenvectors.T
    assert np.allclose((abs(states)**2).sum(axis=1), 1)
    assert np.allclose(np.einsum("ti,ij,tj->t", states.conj(), hamiltonian, states),
                       vacuum @ hamiltonian @ vacuum)
    benchmark = abs(states)**2 @ observables.T
    write_csv("exact_dynamics.csv", ["wt", "density", "flux_1", "flux_2", "flux_3"],
              zip(times, *benchmark.T))
    summary = []
    for noise in (0, 0.001, 0.01):
        for steps in (2, 4, 8, 16, 32):
            predictions, infidelities = [], []
            for time, state in zip(times[1:], states[1:]):
                rho, depth, count = circuit(gates, vacuum, time, steps, noise)
                predictions.append(observables @ rho.diagonal().real)
                infidelities.append(1 - np.vdot(state, rho @ state).real)
            errors = np.abs(np.array(predictions) - benchmark[1:])
            summary.append((noise, steps, depth, count, errors[:, 0].max(),
                            errors[:, 1:].max(), max(infidelities)))
    write_csv("depth_noise.csv", ["noise", "steps", "native_depth", "two_qubit_gates",
              "max_density_error", "max_flux_error", "max_infidelity"], summary)
    peak = benchmark[:, 0].argmax()
    print(f"Peak sampled exact density {benchmark[peak, 0]:.6f} at wt={times[peak]:.2f}")
    for noise, steps, depth, count, density_error, flux_error, infidelity in summary:
        print(f"{noise:5g} {steps:2d} {depth:3d} {count:3d} "
              f"{density_error:.6f} {flux_error:.6f} {infidelity:.6f}")


if __name__ == "__main__":
    main()
