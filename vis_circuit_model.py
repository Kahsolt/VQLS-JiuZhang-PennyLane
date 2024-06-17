#!/usr/bin/env python3
# Author: Armit
# Create Time: 2024/06/17 

# 查看线路对应的数学模型

import numpy as np
import sympy as sp
from sympy import Symbol as S

v0 = np.asarray([1, 0]) 
I = np.asarray([
  [1, 0],
  [0, 1],
])
RY = lambda θ: np.asarray([
  [sp.cos(S(θ, real=True)/2), -sp.sin(S(θ, real=True)/2)],
  [sp.sin(S(θ, real=True)/2),  sp.cos(S(θ, real=True)/2)],
])
CNOT = np.asarray([
  [1, 0, 0, 0],
  [0, 1, 0, 0],
  [0, 0, 0, 1],
  [0, 0, 1, 0],
])
rCNOT = np.asarray([
  [1, 0, 0, 0],
  [0, 0, 0, 1],
  [0, 0, 1, 0],
  [0, 1, 0, 0],
])


if not 'swap-distro':
  '''
  The swap-distro circuit is like:
    --RY--o--RY--x--RY--
          |      |
    ------x------o-----
  项: f*f
  参数: a1/2, (a2+a3)/2, (a2-a3)/2
  '''

  matrix = np.kron(RY('a3'), I) @ rCNOT @ np.kron(RY('a2'), I) @ CNOT @ np.kron(RY('a1'), I)
  state = matrix @ np.kron(v0, v0)

  print('[state]')
  for it in state: print(it)
  print()

  print('[state (simplied)]')
  state = [it.simplify() for it in state]
  for it in state: print(it)
  print()

  print('[state (evalued with pretrained params)]')
  for it in state:
    print(it.evalf(subs={
      S('a1', real=True): 0.782011851189378,
      S('a2', real=True): 0.04758310327698872,
      S('a3', real=True): 0.44237422297673984,
    }))


if 'wtf':
  '''
  The wtf circuit is like:
    --RY--o--RY--
          |
    --RY--x------
  项: f*f*f + g*g*g
  参数: a1/2, a2/2, a3/2
  '''

  matrix = np.kron(RY('a3'), I) @ CNOT @ np.kron(RY('a1'), RY('a2'))
  state = matrix @ np.kron(v0, v0)

  print('[state]')
  for it in state: print(it)
  print()

  print('[state (simplied)]')
  state = [it.simplify() for it in state]
  for it in state: print(it)
  print()

  print('[state (evalued with pretrained params)]')
  for it in state:
    print(it.evalf(subs={
      S('a1', real=True): -0.047208991809220814,
      S('a2', real=True): 0.7808882554867012,
      S('a3', real=True): 0.509390326201012,
    }))
