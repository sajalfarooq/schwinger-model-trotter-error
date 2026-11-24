# Four-site Schwinger model pilot

**Question.** Can a gate-based circuit reproduce the real-time excitation density and electric flux of a four-site Schwinger model? How does the error depend on circuit depth and two-qubit gate noise?

## Method

The pilot uses an open chain, the bare staggered vacuum, zero left boundary flux, and the ratios J/w = 1 and m/w = 0.5 used in Fig. 2(a) of [Martinez et al., Nature (2016)](https://arxiv.org/abs/1605.04570). These ratios define a point of comparison. A quantitative overlay with the published curve remains to be done.

In units with w = 1, the spin Hamiltonian consists of nearest-neighbor hopping w(XX + YY)/2, staggered mass m(-1)^n Z_n/2 summed over one-based sites n, and electric energy J times the sum of squared link fluxes. Gauss's law gives the flux on link n as the boundary flux plus the cumulative staggered charge through site n.

The exact calculation builds the hopping matrix by swapping opposite neighboring bits, then diagonalizes the resulting 16 by 16 Hamiltonian. It checks vacuum charge and flux, total charge conservation, state normalization, and energy conservation. The circuit independently implements hopping as XX and YY rotations and represents the diagonal term using Z and ZZ rotations plus an irrelevant constant. Repeating the ordered rotations gives a first-order product formula.

All-pairs connectivity and native Z, XX, YY, and ZZ rotations define the circuit depth. The noise model applies a two-qubit depolarizing channel after each two-qubit rotation. With probability p, the channel replaces the pair by the maximally mixed state and removes its correlations with the remaining qubits. Density-matrix simulation averages the noise exactly and uses no finite-shot sampling. The data describe simulated circuits under this synthetic noise model.

## Reproduce

With Python 3.10 or newer and NumPy installed, run:

    python schwinger.py

The script writes exact_dynamics.csv at 21 times from wt = 0 to 2 and depth_noise.csv for 15 configurations. Circuit errors are maxima over the 20 nonzero times. Flux error also takes the maximum over the three links.

## Interpretation

At ideal depth 289, the worst sampled absolute errors in excitation density and link flux are approximately 0.00068 and 0.00819, respectively. This gives a concrete accuracy benchmark for digital simulation of a finite lattice gauge theory.

At p = 0.001, excitation-density error is lowest among the tested circuits at depth 73. Link flux instead favors depth 145, so the best depth depends on the observable. At p = 0.01, both favor depth 37 among the tested circuits. The pilot uses one lattice size and parameter point. Results on hardware would additionally depend on compilation, connectivity, calibration, state preparation, and measurement.

## Which approximation buys more accuracy per gate?

A finer time step increases the number of gates exposed to noise. Before assigning a noise probability, we can ask how much physical accuracy each gate buys in the ideal limit. We compared the original sequence with a symmetric sequence that applies half rotations in forward and reverse order. We merged adjacent identical rotations and counted each remaining two-qubit Pauli rotation once. Each method used the largest tested step count that fit a common budget of 144 rotations for four sites and 320 for six. These retrospective comparison points correspond to 16 first-order steps. The full sweep includes five step counts per method.

| Sites | Formula order | Rotations used | Maximum density error | Maximum link flux error |
|---:|---:|---:|---:|---:|
| 4 | 1 | 144 | 0.002739 | 0.017115 |
| 4 | 2 | 136 | 0.007312 | 0.010480 |
| 6 | 1 | 320 | 0.002318 | 0.011086 |
| 6 | 2 | 312 | 0.007127 | 0.006380 |

At both sizes the choice depends on what we want to measure. The symmetric sequence cuts the largest flux error by 39% on four sites and 42% on six. Yet the first-order sequence gives the smaller density error at both budgets. A researcher using only excitation density to choose a circuit would reach a different decision from one using link flux. Comparable counts can have different physical costs on hardware, and the discrete sweep can miss intermediate step counts.

Errors use 81 equally spaced times from wt = 0 to 2. The four-site first-order density error at 16 steps shifts from 0.002730 on the pilot grid to 0.002739, which illustrates the precision limit of a sampled maximum. Both comparisons share mass, coupling, gate order, and open boundaries. The six-site calculation tests sensitivity to lattice length at fixed spacing. A continuum extrapolation would require varying the spacing and controlling the physical volume.

Run python compare_formulas.py after placing it alongside schwinger.py. The script writes all twenty configurations to formula_comparison.csv. The ideal state-vector implementation and the pilot density-matrix implementation agree on the shared four-site density maxima to within 0.0003 across the two time grids. The noisy comparison currently covers first order and four sites only.
