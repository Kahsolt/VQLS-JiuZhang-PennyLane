#!/usr/bin/env python3
# Author: Armit
# Create Time: 2024/06/17 

import numpy as np
from scipy.sparse import csr_matrix
from spinqit import Circuit, Ry, CNOT
from spinqit.algorithm import VQE
from spinqit.algorithm.optimizer import GradientDescent

from utils import preprocess, postprocess, I_, get_fidelity, xv

''' HParam '''
n_qubits = 2
lr = 1.5
iters = 30000
seed = 1024
np.random.seed(seed)

''' Hamiltonian '''
A, b, x = preprocess()
H_A = A.conj().T @ (I_(n_qubits) - b @ b.conj().T) @ A    # Eq. 6
ham = csr_matrix(H_A)

''' Ansatz '''
ansatz = 'swap_distro'
if ansatz == 'two_local':
  param_shape = (6,)
  circ = Circuit()
  qv = circ.allocateQubits(n_qubits)
  p = circ.add_params(shape=param_shape)
  circ << (Ry, qv[0], p[0])
  circ << (Ry, qv[1], p[1])
  circ << (CNOT, [qv[0], qv[1]])
  circ << (Ry, qv[0], p[2])
  circ << (Ry, qv[1], p[3])
  circ << (CNOT, [qv[1], qv[0]])
  circ << (Ry, qv[0], p[4])
  circ << (Ry, qv[1], p[5])
elif ansatz == 'swap_distro':
  param_shape = (3,)
  circ = Circuit()
  qv = circ.allocateQubits(n_qubits)
  p = circ.add_params(shape=param_shape)
  circ << (Ry, qv[0], p[0])
  circ << (CNOT, [qv[0], qv[1]])
  circ << (Ry, qv[0], p[1])
  circ << (CNOT, [qv[1], qv[0]])
  circ << (Ry, qv[0], p[2])

''' Train '''
optim = GradientDescent(maxiter=iters, learning_rate=lr, tolerance=1e-30, verbose=True)
vqe = VQE(ham, optim, ansatz=circ, params=np.zeros(param_shape))
loss_list = vqe.run(mode='spinq', grad_method='param_shift')
losses = [it.item() for it in loss_list]
print('optimized params:', vqe.optimized_params)

''' PMeasure '''
x_tilde = np.asarray([it.real for it in vqe.optimized_result.states])
print(r'|\tilde{x}>:', x_tilde)
print('fid:', get_fidelity(x_tilde, x))
xv_hat = postprocess(x_tilde)
print('x:', xv_hat)
print('L1 err:', np.abs(xv.flatten() - xv_hat).mean())
print()
