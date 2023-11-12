import numpy as np
import random as rnd
import matplotlib.pyplot as plt

def make_state(N):
  return np.ones((N, N))

def H(J, X):
  N = len(X)
  sum = 0.0
  for i in range(N):
    for t in range(N):
      sum += X[i,t] * X[i, (t + 1) % N]
  return -J * sum

def M(X):
  return np.sum(X)

def T(X):
  N = len(X)
  return np.roll(X, 1, axis=1)

def evolve(X, J, T):
  return X * np.exp(-2.0 * J * T)

def run1(N, J, T, nsteps):
  X = make_state(N)
  Hs = []
  Ms = []
  for step in range(nsteps):
    Hs.append(H(J, X))
    Ms.append(M(X))
    X = evolve(X, J, T)
  return Hs, Ms

def plot(Hs, Ms, N, J, T, nsteps):
  f, ax = plt.subplots()
  e_avg = np.average(Hs[:-1000])
  m_avg = np.average(Ms[:-1000])
  ax.plot(Hs, label='H(J, X)')
  ax.plot(Ms, label='M(X)')
  ax.axhline(e_avg, color='r', linestyle='--', label=f"avg={e_avg}")
  ax.axhline(m_avg, color='g', linestyle='--', label= f"avg={m_avg}")
  ax.title.set_text(f"N={N}, J={J}, beta={1/T}, nsteps={nsteps}")
  # ax.text(1.0, 0.0, f"N={N}, J={J}, T={T}, nsteps={nsteps}",
  #    horizontalalignment='right',
  #    verticalalignment='bottom',
  #    transform = ax.transAxes)

  ax.legend()
  plt.show()


def delta_E_flip(J, spin_config, x, y):
  """ Compute the chamge in energy due to a spin flip at (x,y) in a 2D lattice."""
  N = len(spin_config)
  # J = -J * -1 - explicitly showing the flip the flip of spin_config[x,y] * -J
  return J * spin_config[x, y] * (spin_config[(x-1)%N, y] +
                               spin_config[(x+1)%N, y] +
                               spin_config[x, (y-1)%N] +
                               spin_config[x, (y+1)%N])


def run2(N, J, T, nsteps):
  beta = 1.0 / T

  spin_config = make_state(N)

  Es = [H(J, spin_config)]
  print(f"Initial Energy:{Es[0]}")
  Ms = [M(spin_config)]
  print(f"Initial Magnetization:{Ms[0]}")

  cur_E = Es[0]
  for _ in range(nsteps):
    # Flip a single spin chosen at random:
    xndx = rnd.randint(0, N-1)
    yndx = rnd.randint(0, N-1)
  
    new_E = cur_E + delta_E_flip(J, spin_config, xndx, yndx)
    accept = new_E < cur_E or rnd.random() < np.exp(-beta * (new_E - cur_E))
    # if _ < 40: print(f"[{_}] ({xndx},{yndx}) cur_E:{cur_E} new_E:{new_E}  accept:{accept}")
    if accept:
      spin_config[xndx, yndx] *= -1 # flip the spin @(x,y)
      Es.append(new_E)
      Ms.append(M(spin_config))
      cur_E = new_E
    else:
      Es.append(cur_E)
      Ms.append(Ms[-1])

  return Es, Ms

def run1_setup():
  N = 10
  J = 1.0
  T = 1.0
  nsteps = 1000
  return N, J, T, nsteps

def do_run1():
  Hs, Ms = run1(*run1_setup())
  plot(Hs, Ms, *run1_setup())

def run2_setup():
  N=32
  J = 1.0
  T = 2.0
  nsteps = N*N*100
  return N, J, T, nsteps

def do_run2(): 
  Hs, Ms = run2(*run2_setup())
  plot(Hs, Ms, *run2_setup())

if __name__ == '__main__':
  do_run2()
