"""Exact rational coefficients for DLMF 8.12.9-11, generated with Python 3.12.

Lagrange inversion of eta^2/2 = t-log(1+t) gives t(eta). The identity
t'(eta)=1+eta+eta*c0(eta) then gives c0 without subtracting singular terms.
These are defining mathematical coefficients, not fitted probability data.
"""
from fractions import Fraction as F

def multiply(a, b, degree):
    c = [F(0)] * (degree + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:degree + 1 - i]):
            c[i+j] += x*y
    return c

def inverse_coefficient(n):
    degree = n-1
    # (eta/t)^2 = 2*(t-log(1+t))/t^2.
    u = [F(0)] + [F(2*(-1)**j, j+2) for j in range(1, degree+1)]
    power = [F(1)] + [F(0)]*degree
    factor = F(1)
    answer = F(0)
    for m in range(degree+1):
        answer += factor*power[degree]
        power = multiply(power, u, degree)
        factor *= (F(-n, 2)-m)/(m+1)
    return answer/n

c0 = [(n+2)*inverse_coefficient(n+2) for n in range(24)]
c0[0] -= 1
rows = [c0]
for k, g in enumerate([F(1,12), F(1,288), F(-139,51840)], 1):
    previous = rows[-1]
    rows.append([(n+2)*previous[n+2]+(-1)**k*g*c0[n]
                 for n in range(len(previous)-2)])
for row in rows:
    print("new double[] { " + ", ".join(format(float(v), ".17g") for v in row) + " },")
