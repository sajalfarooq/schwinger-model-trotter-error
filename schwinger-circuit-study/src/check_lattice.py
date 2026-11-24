"""Check all four observables with direct Taylor evolution in the neutral sector."""
from itertools import product
from pathlib import Path
import numpy as np

states = [s for s in product([1, -1], repeat=4) if sum(s) == 0]
stagger = np.array([-1, 1, -1, 1])
h = np.zeros((6, 6))
observables = []
for i, s in enumerate(states):
    flux = np.cumsum((s + stagger) / 2)[:-1]
    observables.append([(1 + s * stagger).mean() / 2, *flux])
    h[i, i] = np.dot(s, stagger) / 4 + flux @ flux
    for j in range(3):
        if s[j] != s[j+1]:
            swapped = list(s)
            swapped[j:j+2] = swapped[j:j+2][::-1]
            h[states.index(tuple(swapped)), i] = 1
vac = np.eye(6)[states.index((1,-1,1,-1))]
density = np.diag(np.array(observables)[:,0])
comm = h @ density - density @ h
curvature = -vac @ (h @ comm - comm @ h) @ vac
assert curvature == 3
data = np.loadtxt(Path(__file__).resolve().parents[1] / "data/exact_dynamics.csv",
                  delimiter=",", skiprows=1)
errors = []
for t, *reference in data:
    psi = vac.astype(complex)
    term = psi.copy()
    for order in range(1, 100):
        term = (-1j*t/order)*(h @ term)
        psi += term
        if np.linalg.norm(term) < 1e-17:
            break
    else:
        raise ArithmeticError("Taylor series failed to converge")
    errors.extend(abs(abs(psi)**2 @ observables - reference))
assert max(errors) < 1e-10
print(f"Vacuum curvature {curvature:g}; {len(errors)} comparisons; max error {max(errors):.3g}")
