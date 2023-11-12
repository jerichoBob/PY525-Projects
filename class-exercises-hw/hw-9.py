import numpy as np
import matplotlib as plt
import lib.operator as op
import lib.potential as pt
import lib.riccati as rc
import lib.system as sy


print("="*80)
def gauleg(n, a, b):
  xs, ws = np.polynomial.legendre.leggauss(n)

  xps = np.copy(xs)
  wps = np.copy(ws)

  for i in range(n):
    xps[i] = 0.5 * (xs[i] + 1.0) * (b - a) + a
    wps[i] = 0.5 * ws[i] * (b - a)

  return xps, wps
def PV(a, b, n, sing):
  R1 = 0
  Integral = 0
  R2 = np.log(np.absolute((b-sing)/(a-sing)))
  R3 = np.pi*1j

  nleft = int(n * np.abs(sing-a) / (b-a))
  nright = int(n * np.abs(b-sing) /(b-a))
  xps1, wps1 = gauleg(nleft, a, sing)
  xps2, wps2 = gauleg(nright, sing, b)

  x = np.append(xps1, xps2)
  w = np.append(wps1, wps2)

  for i in range(n-1):
    R1 += w[i]/(sing - x[i])

  R = R1+R2+R3

  Singularity = R*np.exp(-sing)

  for i in range(n-1):
    Integral += w[i] * (np.exp(-x[i])/(x[i]-1))

  Sol = Singularity + Integral
  return Sol


a = 0
b = 3
n = 100
x0 = 1
I = PV(a, b, n, x0)
print(np.real(I))

print("="*80)
def gauleg_pv(n, a, b, x0):
  n1 = max(2, int(n * abs(x0 - a) / abs(b - a)))
  n2 = n - n1

  xs1, ws1 = gauleg(n1, a, x0)
  xs2, ws2 = gauleg(n2, x0, b)

  xs = np.append(xs1, xs2)
  ws = np.append(ws1, ws2)

  R = np.sum([xw[1] / (x0 - xw[0]) for xw in zip(xs, ws)]) + np.log(abs((b - x0) / (a - x0)))

  return np.append(xs, x0), np.append(ws, R)
def pv_integrate(a, b, n, sing):
  R1 = 0
  Integral = 0
  R2 = np.log(np.absolute((b-sing)/(a-sing)))
  R3 = np.pi*1j

  xs, ws = gauleg_pv(n, a, b, sing)

  for i in range(n-1):
    R1 += ws[i]/(sing - xs[i])

  R = R1+R2+R3

  Singularity = R*np.exp(-sing)

  for i in range(n-1):
    Integral += ws[i] * (np.exp(-xs[i])/(xs[i]-1))

  Sol = Singularity + Integral
  return Sol

I = pv_integrate(a, b, n, x0)
print(f"I(realpart): {np.real(I)}")

print("="*80)
sys = sy.System()

V = pt.V_Gauss(sys, -4.0, 2.0)
print(f"V:{V}")

print("="*80)
l = 0
p_i = 0.12
p_j = 0.34
print(f"V.get({l}, {p_i}, {p_j}): {V.get(l, p_i, p_j)}")
print("="*80)

n = 64
Lambda = 16
k = 1
def make_mesh(n, Lambda, k):
  xs, ws = gauleg_pv(n, 0, Lambda, k)
  ws = ws.astype(complex)
  ws[-1] = ws[-1] + (1j*np.pi)
  return xs, ws

def V_pk(n, Lambda, k): #I opened up the library files to make this guess, I literally have no clue if this is right
  xs, ws = make_mesh(n, Lambda, k)
  V_out = []
  for i in range(n+1): #minus one because the last point is the k value singularity
    V_out = np.append(V_out, V.get(l, xs[i], k))
  return V_out


G0 = op.G_0(sys)

k = 1.0
print("residue: ", G0.residue(k))

E = sys.e_from_k(k)
q = 1.1
print(f"G0({E}, {q}): {G0(E, q)}")


def K_ij(n, Lambda, k):
  xs, ws = make_mesh(n, Lambda, k)
  K_out = np.ndarray(shape=(n+1, n+1), dtype= np.cdouble)
  E = sys.e_from_k(k)
  l= 0
  for i in range(n+1):
    for b in range(n+1):
      if b == n or i == n:
        K_out[i,b] = (ws[b]*(xs[b]**2))/(2*np.pi**2) * V.get(l, xs[i], xs[b]) * (G0.residue(k)*1j)
      else:
        K_out[i,b] = (ws[b]*(xs[b]**2))/(2*np.pi**2) * V.get(l, xs[i], xs[b]) * (G0(E, xs[b]))
  return K_out
