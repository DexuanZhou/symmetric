基函数形式：
======
以N=2，d=3为例
****

![](https://latex.codecogs.com/svg.latex?\\left(T_{k_{11}}(x_{11})T_{k_{12}}(x_{12})T_{k_{13}}(x_{13})+T_{k_{11}}(x_{21})T_{k_{12}}(x_{22})T_{k_{13}}(x_{23})\right)\cdot%20\left(T_{k_{21}}(x_{11})T_{k_{22}}(x_{12})T_{k_{23}}(x_{13})+T_{k_{21}}(x_{21})T_{k_{22}}(x_{22})T_{k_{23}}(x_{23})\right)\cdot%20\prod_{i=1}^{N}\prod_{j=1}^d(range^2-x_{ij}^2))

函数解释：
=======
|  函数   | 解释 |
|  :----  | ----  |
| basis_function  | 输入x,k返回基函数的值 |
| d2chebyshev  | 输入x,k返回切比雪夫多项式二阶导函数值 |
| laplace_psi  | 输入x,k返回laplace psi |
| V_psi  | 输入x,k返回Vpsi |
| H_psi  | 输入x,k返回Hpsi |
| nsumk  | 输入N,M返回将M分为N个数的所有满足条件的数组，为找k做准备|
| find_k  | 输入d,N,M返回所有满足条件的k |
| random_choose  | 随机选点，并计算所有psi(x)，以及Hpsi(x) |
| least_squares  | 最小二乘|
| Eigenvalue  | 梯度下降|
| reilaygh_quotient_alpha  | 输入系数和选点，归一化系数，以及计算reilaygh_quotient|
| steep  | 黄金分割求步长|
