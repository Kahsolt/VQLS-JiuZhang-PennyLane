# VQLS-JiuZhang

    Contest solution for 2024第三届“量旋杯”大湾区量子计算挑战营

----

Contest page: [https://quantum-challenge.spinq.cn/competitionDetail/profession](https://quantum-challenge.spinq.cn/competitionDetail/profession)  
Team Name: 啊这是什么吃一口  


### Problem

```
[原题 - 出自《九章算术·方程·8》]
今有卖牛二、羊五，以买一十三豕，有余钱一千；
卖牛三、豕三，以买九羊，钱适足；
卖六羊，八豕，以买五牛，钱不足六百。
问：牛、羊、豕价各几何？
答曰：牛价一千二百，羊价五百，豕价三百。

[线性方程组]
 2 * x + 5 * y - 13 * z = 1000
 3 * x - 9 * y +  3 * z = 0
-5 * x + 6 * y +  8 * z = -600
解：牛 x = 1200, 羊 y = 500, 豕 z = 300
```

赛题本质为给定的**线性方程组求解**，其同解方程组为：

$$
\begin{array}{ll}
\begin{bmatrix}
   2 &  5 & -13 \\
   3 & -9 &   3 \\
  -5 &  6 &   8 \\
\end{bmatrix} \begin{bmatrix}
  12 \\
  5 \\
  3 \\
\end{bmatrix} = \begin{bmatrix}
  10 \\
  0 \\
  -6 \\
\end{bmatrix}
\end{array}
$$

已知的 **量子线性求解器 Quantum Linear-system Solver** 算法流派有：

- HHL / QPE-based
- Adiabatic-based: 绝热演化
- VQLS: 变分线路 ⭐
- Grover-based (?)

考虑到赛题对所用量子门和线路深度的限制，**VQLS** 方法应该是唯一正解 🤔


### VQLS in brief

![VQSL](img/VQLS.png)

- 使用线路 $ U $ 来制备态 $ \left| b \right> $，即 $ \left| b \right> = U \left| 0 \right> $
- 使用可学习的含参线路 $ V(\alpha) $ 来制备态 $ \left| x \right> $，即 $ \left| x(\alpha) \right> = V(\alpha) \left| 0 \right> $
- 将系数矩阵 $ A $ 编码为一个 Oracle (其实是 Block-encoding)，即 $ F(A) $ 以实现矩阵乘法 $ \left| \psi \right> = A \left| x(\alpha) \right> $
- 引入一个辅助比特以实现 Hadamard Test，比较两个态 $ \left| \psi \right>  $ 和 $ \left| b \right> $ 的相似度
- 直觉上需要最小化损失 $ C_G = 1 - ||\left< \psi | b \right>||^2 $，但引入惩罚项来避免 barren plateaus
  - 使用 $ C_L = 1 - || \left< x | H_L | x \right> ||^2 $
  - 其中 $ H_L = A^\dagger U (\mathbb{1} - \frac{1}{n} \sum\limits^n_{j=1} \left| 0_j \right> \left< 0_j \right| \otimes \mathbb{1}_{\bar j}) U^\dagger A $


#### refenrence

- essay & tutorial
  - (2018) Solving Linear Systems of Equations by Using the Concept of Grover's Search Algorithm: An IBM Quantum Experience: [https://arxiv.org/abs/1801.00778](https://arxiv.org/abs/1801.00778)
  - (2019) Variational Quantum Linear Solver: [https://arxiv.org/pdf/1909.05820v4.pdf](https://arxiv.org/pdf/1909.05820v4.pdf)
  - Qiskit VQLS tutorial: [https://github.com/qiskit-community/qiskit-textbook/blob/main/content/ch-paper-implementations/vqls.ipynb](https://github.com/qiskit-community/qiskit-textbook/blob/main/content/ch-paper-implementations/vqls.ipynb)
  - VQLS 变分量子算法解线性方程组: [https://blog.csdn.net/qq_43550173/article/details/121591659](https://blog.csdn.net/qq_43550173/article/details/121591659)
  - VQLS 的 MindQuantum 复现: [https://www.cnblogs.com/liniganma/p/17323717.html](https://www.cnblogs.com/liniganma/p/17323717.html)
  - PaddlePaddle-Quantum VQLS: [https://github.com/PaddlePaddle/Quantum/blob/master/applications/linear_solver/introduction_cn.ipynb](https://github.com/PaddlePaddle/Quantum/blob/master/applications/linear_solver/introduction_cn.ipynb)
  - Adiabatic-Linear-Solver-QPanda: [https://github.com/Kahsolt/Adiabatic-Linear-Solver-QPanda](https://github.com/Kahsolt/Adiabatic-Linear-Solver-QPanda)
  - Solving Equations with Grover's Algorithm: [https://www.iap.uni-jena.de/iapmedia/2321/eqt-lecture3](https://www.iap.uni-jena.de/iapmedia/2321/eqt-lecture3)
- dev framework
  - SpinQit: [https://github.com/SpinQTech/SpinQit](https://github.com/SpinQTech/SpinQit)
    - doc: [https://doc.spinq.cn/doc/spinqit/index.html](https://doc.spinq.cn/doc/spinqit/index.html)
  - Mindquantum: [https://www.mindspore.cn/mindquantum/docs/zh-CN/master/index.html](https://www.mindspore.cn/mindquantum/docs/zh-CN/master/index.html)

----
by Armit
2024/6/13
