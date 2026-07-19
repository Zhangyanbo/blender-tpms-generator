# Blender Gyroid 插件更新规格：从伪四边片改为 48 张真正四边参数片

请直接修改当前 Blender 插件的几何生成核心。插件现有的安装方式、注册方式、面板布局和项目结构可以自行适配；本规格只规定新的数学、网格拓扑和验收标准。

最终目标不是“把三角形索引写成四边形”，而是：

1. 一个常规 Gyroid 立方晶胞由 **48 张真正四边参数片**平铺；
2. 每张参数片的几何轮廓确实有四条边，而不是在三角轮廓上人为增加一个 (180^\circ) 角；
3. 每张参数片都有正方形参数域 ([0,1]^2)；
4. 原来两个三角基本域之间的圆弧公共边成为片内曲线，并从最终网格拓扑中消失；
5. 所有输出面均为四边面，不使用 Marching Cubes、Surface Nets、三角化或体素重网格化；
6. 周期边界严格配对；输出不存在重复面、退化面、翻折面和非流形边。

> **2.5 实现修订：**后续推导得到了更直接的解析构造。每对基本三角域先经公共圆
> 反演展开为真实四边复参数域，再由其四条解析边界直接建立 Coons 方形映射。
> 这一构造取代本文第 5–6 节的初始分片坐标与调和重参数化；运行时不再需要
> `solver_resolution` 或 `quadrature_order` 用户参数。旧节保留为历史验收背景。

---

## 1. 必须删除的旧构造

旧版本把一个曲边三角基本域 (O-P-R) 的圆弧 (P-R) 在中点 (Q) 分成两段，再把 (O,P,Q,R) 当作四个角，使用 Coons 映射从正方形映射过去。

这个构造虽然产生四边形索引，但 (Q) 位于一条光滑圆弧内部，几何内角是 (180^\circ)。因此该“面片”的实际轮廓仍是三角形。必须删除这种伪四边参数片。

新版本不得保留以下几何生成路线：

- 三维标量场等值面提取；
- Marching Cubes / Marching Tetrahedra；
- Surface Nets；
- 把 (Q) 当作第四个几何角的单三角域 Coons 映射；
- 先生成三角网格再做 quad remesh；
- 输出两个表面重叠的周期副本。

---

## 2. 精确 Gyroid 基本三角片

### 2.1 Weierstrass 积分

令

\[
K=\operatorname{EllipticK}(1/4),\qquad
K'=\operatorname{EllipticK}(3/4),
\]

\[
\theta=\arctan(K/K'),
\qquad
\kappa=a\frac{\sqrt{K^2+K'^2}}{KK'},
\]

其中 (a) 是论文中 12 片 bounding cell 的边长；内部计算可取 (a=1)。

定义

\[
I_p(\omega)=
\kappa\int_0^\omega
\frac{\tau^p\,d\tau}
{\sqrt{\tau^8-14\tau^4+1}},
\qquad p=0,1,2.
\]

精确 Gyroid 基本片为

\[
X(\omega)=
\begin{pmatrix}
\Re\left[e^{i\theta}(I_0-I_2)\right]\\[2mm]
\Re\left[e^{i\theta}i(I_0+I_2)\right]\\[2mm]
\Re\left[e^{i\theta}2I_1\right]
\end{pmatrix}.
\]

不要用常见的三角函数等值面

\[
\sin x\cos y+\sin y\cos z+\sin z\cos x=0
\]

代替它。后者只是近似隐式面，不是这里需要的参数曲面。

### 2.2 复平面中的曲边三角域

积分域由正实轴、正虚轴和圆

\[
\left|\omega+\frac{1+i}{\sqrt2}\right|=\sqrt2
\]

围成。

从原点发出的射线

\[
\omega=\rho e^{i\phi},\qquad 0\leq\phi\leq\frac\pi2
\]

与圆弧的交点半径为

\[
\rho_{\max}(\phi)=
\frac{-\sqrt2(\cos\phi+\sin\phi)
+\sqrt{2(\cos\phi+\sin\phi)^2+4}}{2}.
\]

因此可用三角坐标

\[
\omega(r,\phi)=r\,\rho_{\max}(\phi)e^{i\phi},
\qquad 0\le r\le1,
\quad 0\le\phi\le\frac\pi2.
\]

三个几何顶点为

\[
O=0,
\quad
P=\rho_{\max}(0),
\quad
R=i\rho_{\max}(\pi/2).
\]

### 2.3 数值积分

沿径向路径令

\[
\tau=t\omega,\qquad 0\le t\le1,
\]

则

\[
I_p(\omega)=
\kappa\omega^{p+1}
\int_0^1
\frac{t^p\,dt}
{\sqrt{(t\omega)^8-14(t\omega)^4+1}}.
\]

在 (P,R) 附近存在平方根端点奇性。使用

\[
t=2s-s^2,qquad dt=2(1-s)ds,
\]

再做 Gauss–Legendre 求积。建议阶数 160–240。沿每条积分路径对根号内部复数的相位做 `unwrap`，从 (t=0) 的正根连续选择平方根分支。

伪代码：

```python
nodes, weights = leggauss(order)
s = 0.5 * (nodes + 1)
w = 0.5 * weights
t = 2*s - s*s
w = w * 2*(1-s)

tau = t * omega
z = tau**8 - 14*tau**4 + 1
sqrt_z = sqrt(abs(z)) * exp(0.5j * unwrap(angle(z)))

Ip = kappa * omega**(p+1) * sum(w * t**p / sqrt_z)
```

---

## 3. 生成常规晶胞中的 96 个基本三角片

以下矩阵采用行向量实现时，应计算 `point @ A.T + b`。

令

```python
q = 1 / sqrt(2)
h = 1 / 4
```

12 个基本片变换为：

```python
PATCH_TRANSFORMS = [
    ([[-q,-q,0],[-q, q,0],[0,0, 1]], [0,   2*h, 3*h]),
    ([[ q,-q,0],[-q,-q,0],[0,0,-1]], [0,   2*h, 3*h]),
    ([[-q, q,0],[0,0, 1],[ q, q,0]], [2*h, 3*h, 4*h]),
    ([[-q,-q,0],[0,0,-1],[-q, q,0]], [2*h, 3*h, 4*h]),
    ([[0,0, 1],[ q, q,0],[ q,-q,0]], [3*h, 4*h, 2*h]),
    ([[0,0,-1],[-q, q,0],[ q, q,0]], [3*h, 4*h, 2*h]),
    ([[ q, q,0],[ q,-q,0],[0,0,-1]], [4*h, 2*h,   h]),
    ([[-q, q,0],[ q, q,0],[0,0, 1]], [4*h, 2*h,   h]),
    ([[ q,-q,0],[0,0,-1],[-q,-q,0]], [2*h,   h,   0]),
    ([[ q, q,0],[0,0, 1],[ q,-q,0]], [2*h,   h,   0]),
    ([[0,0,-1],[-q,-q,0],[-q, q,0]], [  h,   0, 2*h]),
    ([[0,0, 1],[ q,-q,0],[-q,-q,0]], [  h,   0, 2*h]),
]
```

8 个 bounding-cell 变换为：

```python
BLOCK_TRANSFORMS = [
    ([[ 1,0,0],[0, 1,0],[0,0,1]], [0,0, 0]),
    ([[-1,0,0],[0,-1,0],[0,0,1]], [2,1, 0]),
    ([[-1,0,0],[0, 1,0],[0,0,1]], [1,1, 0]),
    ([[ 1,0,0],[0,-1,0],[0,0,1]], [1,2, 0]),
    ([[ 1,0,0],[0,-1,0],[0,0,1]], [0,1,-1]),
    ([[-1,0,0],[0, 1,0],[0,0,1]], [2,0,-1]),
    ([[-1,0,0],[0,-1,0],[0,0,1]], [1,2,-1]),
    ([[ 1,0,0],[0, 1,0],[0,0,1]], [1,1,-1]),
]
```

对每个基本点 (x)，依次计算

```python
x1 = x @ patch_A.T + patch_b
x2 = x1 @ block_A.T + block_b
```

论文坐标中的完整晶胞为

\[
x''\in[0,2],\qquad y''\in[0,2],\qquad z''\in[-1,1].
\]

最终归一化到 Blender 单位晶胞：

```python
X = x2.x / 2
Y = x2.y / 2
Z = (x2.z + 1) / 2
```

必须按上述名义晶格坐标归一化。不要按照曲面自身的 `min/max bounding box` 分别缩放三个方向，否则会破坏周期平移关系。

---

## 4. 把 96 个三角片天然配成 48 个四边片

### 4.1 配对的几何原则

每个基本三角片有三条边：

1. (O-P) 径向边；
2. (O-R) 径向边；
3. (P-R) 圆弧边。

在完整晶胞中，每条圆弧边只与另一张基本三角片的圆弧边重合。96 条圆弧边因此形成 48 个互不相交的配对：

\[
Q_k=T_{a_k}\cup_{P-R}T_{b_k},
\qquad k=1,\ldots,48.
\]

合并后，公共圆弧 (P-R) 进入 (Q_k) 内部；两个三角片的四条径向边构成 (Q_k) 的四条外边界。因此 (Q_k) 是真正的曲边四边形。

### 4.2 不要硬编码配对编号

从每张三角片的圆弧上等距采样，例如 13–25 个点。对变换后的三维坐标执行：

1. 对坐标取模 1，以识别周期平移等价的曲线；
2. 把接近 1 的坐标归为 0；
3. 坐标舍入到 (10^{-5})；
4. 分别生成正向和反向点序列 key，取字典序较小者作为 canonical key。

```python
def periodic_curve_key(curve):
    c = mod(curve, 1.0)
    c[isclose(c, 1.0, atol=2e-5)] = 0.0
    forward = tuple(round(c, 5).ravel())
    backward = tuple(round(c[::-1], 5).ravel())
    return min(forward, backward)
```

按 key 分组后必须得到：

```text
分组数 = 48
每组曲线数 = 2
```

在这个具体构造中，每一对圆弧的方向相反。配对的两个三角片分别记为 `lower` 与 `upper`。

---

## 5. 从正方形到合并四边域的初始精确映射

令正方形坐标为 ((s,t)\in[0,1]^2)，以对角线 (t=s) 暂时把正方形分为两半。

### 5.1 下半部 (t\le s)

映射到第一张三角片，其三个重心坐标为

\[
\lambda_O=s-t,
\qquad
\lambda_P=1-s,
\qquad
\lambda_R=t.
\]

### 5.2 上半部 (t\ge s)

映射到第二张三角片。由于公共圆弧方向相反，使用

\[
\lambda_O=t-s,
\qquad
\lambda_P=s,
\qquad
\lambda_R=1-t.
\]

### 5.3 重心坐标到复参数

令

\[
r=\lambda_P+\lambda_R,
\qquad
\phi=\frac\pi2
\frac{\lambda_R}{\lambda_P+\lambda_R}.
\]

当 (r=0) 时直接取 (omega=0)；否则

\[
\omega=
r\rho_{\max}(\phi)e^{i\phi}.
\]

分别通过精确 (X(\omega)) 以及两张三角片各自的刚体变换得到三维位置。记这个分片精确映射为

\[
P_0(s,t).
\]

在对角线上，两边位置误差应接近数值积分误差，测试值应小于 (10^{-9})。

但是不要把 (P_0) 直接作为最终参数化。它只有 (C^0) 连续；参数线穿过对角线时会出现可见折角。必须执行下一节的调和重参数化。

---

## 6. 调和重参数化：删除隐藏对角线的参数折角

### 6.1 数学定义

在合并后的光滑曲面 (Q_k) 上求调和坐标

\[
\Delta_{Q_k}u=0,
\qquad
\Delta_{Q_k}v=0,
\]

并把四条外边界分别映射到单位正方形的四条边。

设

\[
H:(s,t)\mapsto(u,v)
\]

是从初始分片坐标到调和方形坐标的映射，则最终四边参数片为

\[
\boxed{
\Phi_k(u,v)=P_{0,k}\!\left(H^{-1}(u,v)\right)
},
\qquad (u,v)\in[0,1]^2.
\]

所有 48 个 (Q_k) 彼此刚体等距，因此只需要计算一次 (H^{-1})，其余 47 片复用同一重参数化。

### 6.2 离散实现

建议 `solver_resolution = 44`，可在 32–64 之间设置并缓存结果。

1. 在初始 ((s,t)) 正方形上建立 ((N+1)\times(N+1)) 规则采样；
2. 用 (P_0(s,t)) 计算一个代表性宏观片的精确三维坐标；
3. 仅为了求解参数坐标，把每个计算单元沿 `(i,j) -> (i+1,j+1)` 分为两个辅助三角形；
4. 计算三维三角网格的 cotangent Laplacian：

\[
w_{ij}=\frac12(\cot\alpha_{ij}+\cot\beta_{ij}),
\]

\[
L_{ii}=\sum_jw_{ij},
\qquad
L_{ij}=-w_{ij};
\]

5. 边界顶点的 ((u,v)) 固定为正方形四边上的对应位置；
6. 对内部顶点分别求解

\[
L_{II}u_I=-L_{IB}u_B,
\qquad
L_{II}v_I=-L_{IB}v_B;
\]

7. 检查调和 UV 中所有辅助三角形的有向面积严格为正，不允许翻转；
8. 在调和 UV 平面中定位最终规则 ((u,v)) 网格点所在的辅助三角形；
9. 用重心插值求出相应的初始 ((s,t)=H^{-1}(u,v))；
10. 再通过精确 (P_0(s,t)) 计算最终三维坐标。

辅助三角形只能存在于调和坐标求解器内部。不得把它们传给 Blender mesh，也不得把对角线作为最终边输出。

### 6.3 缓存

`H^{-1}` 只依赖精确基本片的内蕴度量和求解分辨率，不依赖晶胞尺寸、刚体方向或宏观片编号。插件启动后第一次生成时计算并缓存：

```python
cache_key = (solver_resolution, quadrature_order)
```

当用户只修改最终 `quad_subdivisions`、晶胞缩放或显示参数时，不重新求解调和坐标。

---

## 7. Blender 最终网格生成

### 7.1 参数

插件至少保留或增加以下参数：

```text
quad_subdivisions: 每张宏观片每个方向的四边形数量
    1 -> 48 个四边面，仅表示宏观拓扑，几何非常粗
    2 -> 192 个四边面，建议默认值
    4 -> 768 个四边面，较平滑
    8 -> 3072 个四边面

solver_resolution: 调和参数求解分辨率，默认 44
quadrature_order: Weierstrass 积分阶数，默认 200 或 240
cell_scale: Blender 中晶胞边长
```

### 7.2 顶点和面

对每张宏观片 (k=1,\ldots,48)，取

\[
u_i=i/n,qquad v_j=j/n,
\qquad i,j=0,\ldots,n,
\]

并计算

\[
p_{kij}=\Phi_k(u_i,v_j).
\]

四边面索引固定为

```python
(i, j) -> (i+1, j) -> (i+1, j+1) -> (i, j+1)
```

即：

```python
face = (
    index[i,   j],
    index[i+1, j],
    index[i+1, j+1],
    index[i,   j+1],
)
```

相邻宏观片的公共边顶点必须合并。建议用

```python
key = tuple(round(point / 3e-6).astype(int))
```

做空间哈希。

不要把相差整数晶格向量的周期边界顶点在单个晶胞对象中合并，因为它们位于立方体的相对两侧。它们只在周期商空间中等价；复制晶胞后再按平移位置合并。

### 7.3 法向一致性

构建 edge-to-face 邻接表，对所有面做 BFS 定向：任意一条内部公共边在两个相邻面中的有向顺序必须相反。完成后再根据需要整体反转所有面，以选择 Gyroid 的一侧作为正法向。

### 7.4 Blender API

最终必须直接创建 quad polygon：

```python
mesh.from_pydata(vertices, [], quad_faces)
mesh.update()
```

`quad_faces` 中每项长度必须严格等于 4。不要在几何生成之后调用会自动三角化或重新拓扑的操作。

平滑显示只能改变 shading，不能改变拓扑。例如可以使用 smooth polygon normals，但不要用 Remesh modifier 生成最终几何。

---

## 8. 周期边界和拓扑检查

### 8.1 单个立方晶胞

单个开边界晶胞中：

- 内部边邻接两个四边面；
- 晶胞切口上的边邻接一个四边面；
- 不得存在邻接三个或更多面的边。

### 8.2 周期商空间

验证周期性时，把顶点映射为

```python
periodic_point = mod(point, 1.0)
periodic_point[isclose(periodic_point, 1.0, atol=2e-5)] = 0.0
key = tuple(round(periodic_point, 5))
```

在周期商空间中必须满足：

1. 每条边恰好邻接两个四边面；
2. 所有晶胞边界边都存在唯一的周期配对；
3. Euler 特征

\[
\chi=V-E+F=-8;
\]

4. 对应 genus

\[
g=1-\frac{\chi}{2}=5.
\]

---

## 9. 必须自动执行的验收测试

### 9.1 宏观拓扑

```text
fundamental_triangles = 96
circular_edge_groups = 48
triangles_per_group = 2
genuine_quad_patches = 48
```

### 9.2 不同细分等级

对于 `quad_subdivisions = n`：

\[
F=48n^2.
\]

必须得到：

```text
n = 1 -> 48 quad faces
n = 2 -> 192 quad faces
n = 4 -> 768 quad faces
n = 8 -> 3072 quad faces
```

### 9.3 `n=2` 的参考结果

归一化单位晶胞、空间焊接容差 (3\times10^{-6}) 时，参考结果为：

```text
vertices: 237
quad_faces: 192
non_quad_faces: 0
duplicate_faces: 0
degenerate_quads: 0
folded_quads: 0
nonmanifold_edges: 0
all_periodic_boundary_edges_paired: True
quotient_all_edges_incident_to_two_quads: True
quotient_euler_characteristic: -8
periodic_quotient_genus: 5
```

顶点数允许因焊接实现和浮点容差有极小差别，其余布尔和拓扑指标不允许变化。

### 9.4 几何检查

每个四边面用两条对角三角形仅做验证：

```python
n1 = cross(p1-p0, p2-p0)
n2 = cross(p2-p0, p3-p0)
```

必须满足

```python
dot(n1, n2) > 0
norm(n1) > area_tolerance
norm(n2) > area_tolerance
```

辅助三角形不能进入输出。

### 9.5 隐藏公共边检查

随机选择若干宏观片，检查：

- 两张原三角片在隐藏对角线上的位置误差小于 (10^{-9})；
- 最终 Blender 边集合中不存在该对角线；
- 最终规则参数线穿过原对角线时无可见折角；
- 宏观片外边界恰好为四条连续曲线。

---

## 10. 高阶 CAD 表示

Blender mesh 只是 (Phi_k(u,v)) 的采样。真正的高阶拓扑始终是 48 张方形参数片。

需要更高精度时，对每张调和参数片拟合 tensor-product Chebyshev、Bernstein 或 B-spline：

\[
\Phi_{k,m}(u,v)=
\sum_{a=0}^{m}\sum_{b=0}^{m}
c^{(k)}_{ab}T_a(2u-1)T_b(2v-1).
\]

相邻片共享整条边界控制点；周期边界控制点相差对应晶格向量。提高多项式次数 (m) 可以改善几何精度，而宏观片数量保持为 48。

不要把“提高精度”等同于无限增加拓扑面数。插件中的核心对象应当是

```text
48 x square parameter charts
```

而不是体素等值面或三角网格。

---

## 11. 完成定义

只有同时满足以下条件，更新才算完成：

- 旧的 96 个伪四边三角轮廓被 48 个真实四边宏观片替代；
- 每个宏观片使用调和方形参数化；
- 原圆弧公共边不出现在 Blender 网格边集合中；
- 默认 `quad_subdivisions=2` 时输出 192 个纯四边面；
- 无重复、退化、翻折和非流形结构；
- 周期商空间中每条边邻接两个面，Euler 特征为 (-8)；
- 插件可以在不改变 48 片宏观拓扑的前提下提高采样或高阶拟合精度。
