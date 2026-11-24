"""Independent free-limit check in 1+1 dimensions, metric (+,-)."""
from fractions import Fraction
from math import factorial

import numpy as np

gamma0 = np.array([[1, 0], [0, -1]])
gamma1 = np.array([[0, 1], [-1, 0]])
# Compare every coefficient of (E gamma0 - k gamma1)^2.
assert np.array_equal(gamma0 @ gamma0, np.eye(2))
assert np.array_equal(gamma1 @ gamma1, -np.eye(2))
assert not (gamma0 @ gamma1 + gamma1 @ gamma0).any()
print("Normalized trace of squared slashed momentum: E^2 - k^2")

# Convolve the Taylor coefficients of sin(x) exactly through x^4.
sine = [Fraction(0), Fraction(1), Fraction(0), Fraction(-1, factorial(3)), Fraction(0)]
square = [sum(sine[j] * sine[i-j] for j in range(i+1)) for i in range(5)]
assert square == [0, 0, 1, 0, Fraction(-1, 3)]
print("Squared free band: m^2 + k^2 - a^2 k^4/3 + O(a^4)")
