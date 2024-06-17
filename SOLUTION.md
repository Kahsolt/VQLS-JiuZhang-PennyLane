### Variational Algorithms for Linear Algebra: [arXiv:1909.03898](https://arxiv.org/abs/1909.03898)

> VALA = VQE

- 对任意矩阵 $ M $ 和单位向量 $ \left| v \right> $，求矩阵-向量乘法并规范化 $ \left| v_M \right> = \frac{M \left| v \right>}{||M \left| v \right>||} $
  - 构造哈密顿量 $ H_M = I - \frac{M \left| v \right> \left< v \right| M^\dagger}{||M \left| v \right>||^2} $
  - 则 $ \left| v_M \right> $ 是 $ H_M $ 的基态，且能量为 0
  - 可用 绝热态制备、相位估计 或 **VQE** 等方法去求此基态
- 对于线性方程组求解问题，我们希望计算 $ \left| v_{M^{-1}} \right> = \frac{M^{-1} \left| v \right>}{||M^{-1} \left| v \right>||} $
  - 构造哈密顿量 $ H_{M^{-1}} = M^\dagger (I - \left| v \right> \left< v \right|) M $
  - 则 $ \left| v_{M^{-1}} \right> $ 是 $ H_{M^{-1}} $ 的基态，且能量为 0
  - 注：M 非厄密时可用 [HHL 中的标准构造法](#appendix-a) 扩张为厄密阵，但似乎这个算法并不要求 M 厄密（？
- Ansatz 设计: RZ-RY + cyclic(CNOT + RZ-RY)*n_repeat + rev-cyclic(CNOT)

---

### Variational Quantum Linear Solver: [arXiv:1909.05820](https://arxiv.org/abs/1909.05820v4)

> VQLS = VQE + Hadamard Test

![VQLS](img/VQLS.png)

- 使用线路 $ U $ 来制备参考态 $ \left| b \right> = U \left| 0 \right> $
- 使用含参线路 $ V(\alpha) $ 来制备目标态 $ \left| x \right> = V(\alpha) \left| 0 \right> $
  - Ansatz 设计: centric(CZ, RY)
- 将系数矩阵 $ A $ 编码为一个 Oracle 线路 $ F(A) $ 以实现矩阵乘法 $ \left| \psi \right> \sim A \left| x \right> $
  - 系数矩阵 $ A $ 必须被表达为 LCU，因此必为厄密的 ⚠
  - 由线性性 $ A \left |x \right> = \sum_l A_l \left |x \right> $ 拆分实现，其中 $ A_l $ 是酉的
  - F(A) 是一组可拆装的线路，而非一整个的 PREPARE-SELECT 结构
- 引入一个辅助比特以实现 Hadamard Test，实现态向酉阵的期望投影 $ \left< \psi | U | \psi \right> $
- 损失函数定义得比较复杂，但实际就是 VQE 中求哈密顿量期望
  - 直觉上需要最小化损失 $ C_G = 1 - ||\left< \psi | b \right>||^2 \sim \left< x | H_G | x \right> $, 其中 $ H_G = A^\dagger (1 - \left| b \right> \left< b \right|) A $
  - 也可以使用 "local" 版本 $ C_L = \left< x | H_L | x \right> $ 来避免 barren plateaus，其中 $ H_L = A^\dagger U (\mathbb{1} - \frac{1}{n} \sum\limits^n_{j=1} \left| 0_j \right> \left< 0_j \right| \otimes \mathbb{1}_{\bar j}) U^\dagger A $，即每次只关心 $ \left| b \right> $ 中的某一个分量
- 优化收敛后即有 $ V(\hat \alpha) \left| 0 \right> = \left| \hat x \right> \simeq \left| x \right> $

----

### Variational Quantum Linear Solver with Dynamic Ansatz: [arXiv:2107.08606](https://arxiv.org/abs/2107.08606)

> VQLS-DA = VQLS + Dynamic Ansatz

- 使用 VQLS 中的 "global" 版本损失函数 $ C_G $, 因为简单
- Ansatz 在优化中逐渐增加层数
  - Ansatz 设计: HEA(RY, CNOT, linear) := [RY + linear(CNOT)]*n_repeat

----

#### Appendix A

若原始方程组的系数矩阵 $ A $ 不是厄密的，可按 HHL 中的方法构造增广方程组 $ \tilde A \tilde x = \tilde b $：

$$
\begin{array}{ll}
\begin{bmatrix}
   0 & A \\
   A^\dagger & 0 \\
\end{bmatrix} \begin{bmatrix}
  0 \\
  x \\
\end{bmatrix} = \begin{bmatrix}
  b \\
  0 \\
\end{bmatrix}
\end{array}
$$

其中 $ \tilde A $ 是厄密的，可用于构造哈密顿量等。

#### Appendix B

由于量子线路只能制备出单位化的解向量 $ \left| x \right> = \frac{x}{||x||} $，我们用一些小 trick 来恢复解的模长 $ ||x|| $

⚪ 预处理

扩张原方程，在解向量末尾增加一个 scaling indicator，元素值固定为 1，即

$$
\begin{array}{ll}
\begin{bmatrix}
   A & 0 \\
   0 & c \\
\end{bmatrix} \begin{bmatrix}
  x \\
  1 \\
\end{bmatrix} = \begin{bmatrix}
  b \\
  c \\
\end{bmatrix}
\end{array}
$$

其中 $ c \neq 0 $ 为任意常数

⚪ 后处理

得到解的量子态后，读出的振幅应该具有如下形式

$$
\begin{bmatrix}
  \tilde x \\
  s \\
\end{bmatrix}
$$

则模长恢复的解为 $ \hat x = \tilde x / s $，进一步丢掉 scaling indicator 即得原方程组代数解
