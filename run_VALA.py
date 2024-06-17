#!/usr/bin/env python3
# Author: Armit
# Create Time: 2024/06/17 

# impl [arXiv:1909.03898] Variational Algorithms for Linear Algebra
# https://pennylane.ai/qml/demos/tutorial_vqe/
# https://pennylane.ai/blog/2022/06/how-to-choose-your-optimizer/

from pprint import pprint
from utils import *     # allow shadowing
import pennylane as qml
from pennylane import numpy as np
from pennylane.tape import QuantumTape
from pennylane.operation import Operation
import matplotlib.pyplot as plt


''' HParam '''
circ_type = 'simple'
depth = 1
n_qubits = 2
if circ_type == 'original':
  n_param = n_qubits * (2 + 4 * depth)
elif circ_type == 'simple':
  n_param = n_qubits * (1 + 2 * depth)
n_shots = 10000
lr = 0.8
iters = 500

print('depth:', depth)
print('n_qubits:', n_qubits)
print('n_param:', n_param)
print('lr:', lr)
print('iters:', iters)
print()


''' Data '''
A, b, x = preprocess()
print_matrix(A, 'A')
print_matrix(b, '|b>')
print_matrix(x, '|x>')
print()


''' Hamiltonian '''
nq = int(np.ceil(np.log2(A.shape[0])))
assert n_qubits == nq
H_A = A.conj().T @ (I_(nq) - b @ b.conj().T) @ A    # Eq. 6
print_matrix(H_A, 'H_A')
H = qml.pauli_decompose(H_A)
print(f'[ham] n_terms: {len(H)}')
print(repr(H))
print()


''' Ansatz '''
sim = qml.device('lightning.qubit', wires=nq)
dev = qml.device('lightning.qubit', wires=nq, shots=n_shots)

def circuit_original(param:ndarray) -> List[Operation]:    # Fig. 1
  ops: List[Operation] = []
  pid = 0
  for i in range(n_qubits):
    ops.append(qml.RZ(param[pid], wires=i)) ; pid += 1
    ops.append(qml.RY(param[pid], wires=i)) ; pid += 1
  for d in range(depth):
    for i in range(n_qubits):
      j = (i + 1) % n_qubits
      ops.append(qml.CNOT(wires=[i, j]))
      ops.append(qml.RZ(param[pid], wires=i)) ; pid += 1
      ops.append(qml.RY(param[pid], wires=i)) ; pid += 1
      ops.append(qml.RZ(param[pid], wires=j)) ; pid += 1
      ops.append(qml.RY(param[pid], wires=j)) ; pid += 1
  for i in reversed(range(n_qubits)):
    j = (i + 1) % n_qubits
    ops.append(qml.CNOT(wires=[i, j]))
  return ops

def circuit_simple(param:ndarray) -> List[Operation]:    # just enough for the concrete problem :)
  ops: List[Operation] = []
  pid = 0
  for i in range(n_qubits):
    ops.append(qml.RY(param[pid], wires=i)) ; pid += 1
  for i in range(n_qubits):
    j = (i + 1) % n_qubits
    ops.append(qml.CNOT(wires=[i, j]))
    ops.append(qml.RY(param[pid], wires=i)) ; pid += 1
    ops.append(qml.RY(param[pid], wires=j)) ; pid += 1
  return ops

circuit = eval(f'circuit_{circ_type}')

@qml.qnode(sim)
def circuit_exp(param:ndarray):
  global circuit
  circuit(param)
  return qml.expval(H)

@qml.qnode(sim)
def circuit_state(param:ndarray):
  global circuit
  circuit(param)
  return qml.state()

@qml.qnode(dev)
def circuit_sample(param:ndarray):
  global circuit
  circuit(param)
  return qml.sample()


''' Train '''
p = np.zeros([n_param], requires_grad=True)
param_list = [p]
loss_list = [circuit_exp(p)]
opt = qml.MomentumOptimizer(stepsize=lr)
for i in range(iters):
  p, loss = opt.step_and_cost(circuit_exp, p)
  loss_list.append(circuit_exp(p))
  param_list.append(p)

  if i % 10 == 0:
    print(f'[{i}/{iters}] loss: {loss_list[-1]}')

p_opt = param_list[-1]
print(f'final loss: {loss_list[-1]}')
print(f'final param: {param_list[-1]}')
print()


''' Plot '''
print(qml.draw(circuit_exp)(p_opt))

plt.plot(loss_list, 'b', alpha=0.75, label='loss')
plt.legend()
plt.tight_layout()
fp = LOG_PATH / f'{Path(__file__).stem}.png'
print(f'>> save loss curve to: {fp}')
plt.savefig(fp, dpi=400)
print()


''' Export QASM '''
ops = circuit(p_opt)
meas = [qml.state()]
tape = QuantumTape(ops, meas)
pprint(tape.circuit)
fp = LOG_PATH / f'{Path(__file__).stem}.qasm'
print(f'>> export QASM to: {fp}')
with open(fp, 'w', encoding='utf-8') as fh:
  fh.write(tape.to_openqasm(measure_all=False))
print()


''' PMeasure '''
x_tilde = circuit_state(p_opt).real
print(r'|\tilde{x}>:', x_tilde)
print('fid:', get_fidelity(x_tilde, x))
xv_hat = postprocess(x_tilde)
print('x:', xv_hat)
print('L1 err:', np.abs(xv.flatten() - xv_hat).mean())
print()


''' QMeasure '''
raw_samples = circuit_sample(p_opt)
samples = []
for sam in raw_samples:
  samples.append(int(''.join(str(bs) for bs in sam), base=2))
probs = np.bincount(samples) / n_shots
x_tilde_approx = np.sqrt(probs)
print(r'|\tilde{x_q}>:', x_tilde_approx)
print('fid:', get_fidelity(x_tilde_approx, x))
xv_hat_approx = postprocess(x_tilde_approx)
print('x_q:', xv_hat_approx)
print('L1 err:', np.abs(xv.flatten() - xv_hat_approx).mean())
print()
