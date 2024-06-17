#!/usr/bin/env python3
# Author: Armit
# Create Time: 2024/06/14 

import random
from pathlib import Path
from typing import List, Tuple

import numpy as np
from numpy import ndarray
from numpy.linalg import norm, inv, eig, eigvals

BASE_PATH = Path(__file__).parent
IMG_PATH = BASE_PATH / 'img' ; IMG_PATH.mkdir(exist_ok=True)
LOG_PATH = BASE_PATH / 'log' ; LOG_PATH.mkdir(exist_ok=True)

''' Const '''

# equation
Am = np.asarray([
  [ 2,  5, -13],
  [ 3, -9,   3],
  [-5,  6,   8],
])
bv = np.asarray([[10, 0, -6]]).T
xv = np.asarray([[12, 5, 3]]).T   # classical solution target
assert np.allclose(Am @ xv, bv)

# states & paulis
v0 = np.asarray([[1, 0]]).T   # |0>
v1 = np.asarray([[0, 1]]).T   # |1>
h0 = np.asarray([[1,  1]]).T / np.sqrt(2)   # |+>
h1 = np.asarray([[1, -1]]).T / np.sqrt(2)   # |->
I_ = lambda nq: np.eye(2**nq)
I = np.asarray([   # pauli-i
  [1, 0],
  [0, 1],
])
X = np.asarray([   # pauli-x
  [0, 1],
  [1, 0],
])
Y = np.asarray([   # pauli-y
  [0, -1j],
  [1j, 0],
])
Z = np.asarray([   # pauli-z
  [1, 0],
  [0, -1],
])
H = np.asarray([
  [1,  1],
  [1, -1],
]) / np.sqrt(2)
RY = lambda θ: np.asarray([   # e^(-i*Y*θ/2)
  [np.cos(θ/2), -np.sin(θ/2)],
  [np.sin(θ/2),  np.cos(θ/2)],
])


''' Matrix Utils '''

is_posdef = lambda A: all([ev > 0 for ev in np.linalg.eigvals(A)])
is_hermitian = lambda H: np.allclose(H, H.conj().T)
is_unitary = lambda U: np.allclose((U.conj().T @ U).real, np.eye(U.shape[0]), atol=1e-6)

def assert_hermitian(H:ndarray):
  assert is_hermitian(H), 'matrix should be hermitian'

def assert_unitary(U:ndarray):
  assert is_unitary(U), 'matrix should be unitary'

def spectral_norm(A:ndarray) -> float:
  '''
  spectral norm (p=2) for matrix 
    non-square: ||A||2 = sqrt(λmax(A*A)) sqrt(A*A的最大特征值) = σmax(A) 最大奇异值
    square: ||A||2 = λmax 最大特征值
  '''
  evs = (np.linalg.eigvalsh if is_hermitian(A) else np.linalg.eigvals)(A)
  return max(np.abs(evs))

def condition_number(A:ndarray) -> float:
  '''
  https://www.phys.uconn.edu/~rozman/Courses/m3511_18s/downloads/condnumber.pdf
  NOTE: condition_number(A) = 1 for unitary A
  '''
  evs = (np.linalg.eigvalsh if is_hermitian(A) else np.linalg.eigvals)(A)
  λmax = max(np.abs(evs))
  λmin = min(np.abs(evs))
  return λmax / λmin

def print_matrix(A:ndarray, name:str='A'):
  if len(A.shape) == 2 and min(A.shape) > 1:
    print(f'{name}: (norm={spectral_norm(A):.4g}, κ={condition_number(A):.4g}, shape={A.shape})')
  else:
    print(f'{name}: (norm={np.linalg.norm(A):.4g}, shape={A.shape})')
  print(A.round(4))


''' QState & Bloch Utils '''

Stat = ndarray
Ham = ndarray

def rand_state() -> Stat:
  psi = np.empty([2, 1], dtype=np.complex64)
  psi.real = np.random.uniform(size=psi.shape)
  psi.imag = np.random.uniform(size=psi.shape)
  psi /= np.linalg.norm(psi)
  return psi

def rand_hermitian(N:int) -> Ham:
  A = np.empty([N, N], dtype=np.complex64)
  A.real = np.random.uniform(size=A.shape)
  A.imag = np.random.uniform(size=A.shape)
  H = A.conj().T @ A    # make it hermitian!
  return H

def eigen_state_of_matrix(A:ndarray, id:int=-1) -> ndarray:
  eigvecs = np.linalg.eig(A)[1]
  if id < 0:
    n_eigvecs = eigvecs.shape[-1]
    id = random.randrange(n_eigvecs)
  return np.expand_dims(eigvecs[:, id], -1)   # [N, 1]

def state_str(psi:Stat) -> str:
  psi = drop_gphase(psi)
  a = psi[0].item().real
  c = psi[1].item().real
  d = psi[1].item().imag
  sign1 = '+' if c >= 0 else '-'
  if sign1 == '-': c, d = -c, -d
  sign2 = '+' if d >= 0 else '-'
  if sign2 == '-': d = -d
  return f'{a:.3f} |0> {sign1} ({c:.3f} {sign2} {d:.3f}i) |1>'

def state_vec(psi:Stat) -> List[complex]:
  return psi.T[0].round(4).tolist()

def drop_gphase(psi:Stat) -> Stat:
  return psi * (psi[0].conj() / np.abs(psi[0]))

def amp_to_bloch(psi:Stat) -> Tuple[float, float]:
  psi = drop_gphase(psi).T[0]
  tht = np.arccos(psi[0].real)
  phi = np.angle(psi[1])
  return tht.item(), phi.item()

def state_norm(psi:ndarray) -> Stat:
  return psi / np.linalg.norm(psi)

def get_fidelity(psi:Stat, phi:Stat) -> float:
  return np.abs(np.dot(psi.conj().T, phi)).item()


''' Equantion '''

def preprocess(hermitize_A:bool=False) -> Tuple[ndarray, ndarray, ndarray]:  # A, b, x
  global Am, bv, xv
  # MAGIC: we add a scaling indicator 1 at the right-bottom of the expanded space
  # it tells us how to rescale the normalized quantum solution back to real :)
  Am_ex = np.asarray([
    [ 2,  5, -13, 0],
    [ 3, -9,   3, 0],
    [-5,  6,   8, 0],
    [ 0,  0,   0, 1],
  ])
  bv_ex = np.asarray([[10, 0, -6, 1]]).T
  xv_ex = np.asarray([[12, 5, 3, 1]]).T
  assert np.allclose(Am_ex @ xv_ex, bv_ex)

  if hermitize_A:
    def make_hermitian(A:ndarray) -> ndarray:
      ''' the general way to make a matrix hermitian in HHL paper '''
      N = A.shape[0]
      Ah = np.zeros([N*2, N*2])
      Ah[:N, N:] = A
      Ah[N:, :N] = A.conj().T
      return Ah

    Am_ex = make_hermitian(Am_ex) ; assert_hermitian(Am_ex)
    bv_ex = np.concatenate([bv_ex, np.zeros_like(bv_ex)])
    xv_ex = np.concatenate([np.zeros_like(xv_ex), xv_ex])

  # normalize
  A = Am_ex / np.linalg.norm(bv_ex)
  b = bv_ex / np.linalg.norm(bv_ex)   # |b>
  x = xv_ex / np.linalg.norm(xv_ex)   # |x> = x / |x|, quantum solution target (up to a GPhase)

  return A, b, x

def postprocess(x_tilde:ndarray) -> ndarray:
  global xv
  # MAGIC: decode the scaling indicator
  x_hat = x_tilde / x_tilde[-1].item()
  x_hat = x_hat[:len(xv)]
  return x_hat


if __name__ == '__main__':
  A, b, x = preprocess()
  print_matrix(A, 'A')
  print_matrix(b, '|b>')
  print_matrix(x, '|x>')

  sn = spectral_norm(A)     # 0.1447793025511878
  κ = condition_number(A)   # 97.22474957798339
  print('sn:', sn)
  print('κ:', κ)

  # test fidelity precision: A|x> -> |b>
  assert np.isclose(get_fidelity(state_norm(A @ x), b), 1.0)
  # test numerical precision: |x> -> x
  assert np.allclose(postprocess(x), xv)

  from code import interact
  interact(local=globals())
