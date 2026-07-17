# TPMS Generator

**Exact** triply periodic minimal surfaces for Blender 4.2+ — Gyroid,
Schwarz P and Schwarz D as clean, all-quad meshes from their analytic
Enneper–Weierstrass parametrizations.

Blender 4.2+ 扩展插件 —— 由解析 Enneper–Weierstrass 参数化生成的**精确**
三周期极小曲面（Gyroid、Schwarz P、Schwarz D），输出干净的全四边形网格。

| Gyroid (2×2×2) | Schwarz P (2×2×2) | Schwarz D (2×2×2) |
|:---:|:---:|:---:|
| ![Gyroid](docs/gyroid_2x2x2.png) | ![Schwarz P](docs/schwarz_p_2x2x2.png) | ![Schwarz D](docs/schwarz_d_2x2x2.png) |

One unit cell, quads following the surface's natural parameter lines —
every vertex lies on the true minimal surface / 单个晶胞，四边形沿曲面天然
参数线分布，每个顶点都严格落在真正的极小曲面上：

| Gyroid cell | Schwarz P cell | Schwarz D cell |
|:---:|:---:|:---:|
| ![Gyroid quad flow](docs/gyroid_quadflow.png) | ![Schwarz P cell](docs/schwarz_p_cell.png) | ![Schwarz D cell](docs/schwarz_d_cell.png) |

[English](#english) | [中文](#中文)

---

## English

The surfaces are evaluated directly from their exact Weierstrass
representations: every vertex lies on the true minimal surface (mean
curvature exactly zero), the quads follow the natural parameter lines,
and the mesh is exactly periodic, so tiled cells weld seamlessly. The
familiar level-set approximation `sin x cos y + sin y cos z +
sin z cos x = 0` appears here only as a validation reference.

### Usage

1. Install the add-on (Edit → Preferences → Add-ons → Install from Disk).
2. 3D viewport N-panel → **TPMS** tab.
3. Choose *Type* (Gyroid / Schwarz P / Schwarz D), set *Cell Size*,
   *Cells X/Y/Z*, *Resolution* → **Generate TPMS**.

The result is one unit-cell mesh plus three Array modifiers. Cell counts
are editable live on the modifiers; apply them when you want a single
mesh. Because the cell boundary is periodic to ~1e-9 of the cell size,
the merged tiles are seamless and watertight across cell boundaries.

- **Resolution** — quads per fundamental-patch edge. One cell =
  `patches × res²` quads with 96 / 48 / 192 patches for Gyroid / P / D
  (res 8 → 6144 / 3072 / 12288). Vertices are exact at any resolution;
  raise res only for smoother silhouettes.
- **Smooth Shading** — uses exact analytic normals (the Gauss map of the
  Weierstrass data).

Note on block edges: the cell is one exact translational unit built from
*whole* fundamental patches. Schwarz D patches
end exactly on the cell faces and the Gyroid's overhang is negligible, so
their blocks look box-clean. Schwarz P patches straddle the cell faces
(P embeds straight lines such as `(t, t+1/2, 1/4)·L` that pierce them), so
a finite P block has a ragged skin even though the tiling itself is
exactly seamless — trim with a Boolean if you need a flat-cut block.

### The mathematics

Based on the exact-computation series by Klinowski and co-workers:

- **Gyroid** — P.J.F. Gandy, J. Klinowski, *Exact computation of the
  triply periodic G ('Gyroid') minimal surface*, Chem. Phys. Lett.
  **321** (2000) 363–371.
  [doi:10.1016/S0009-2614(00)00373-0](https://doi.org/10.1016/S0009-2614(00)00373-0)
  ([PDF](https://mathcurve.com/surfaces.gb/Gyroide/sdarticle%20gyroid.pdf))
- **Schwarz P** — P.J.F. Gandy, J. Klinowski, *Exact computation of the
  triply periodic Schwarz P minimal surface*, Chem. Phys. Lett. **322**
  (2000) 579–586.
  [doi:10.1016/S0009-2614(00)00453-X](https://doi.org/10.1016/S0009-2614(00)00453-X)
  ([PDF](https://www.mathcurve.com/surfaces.gb/schwarz/sdarticle.pdf))
- **Schwarz D** — P.J.F. Gandy, D. Cvijović, A.L. Mackay, J. Klinowski,
  *Exact computation of the triply periodic D ('diamond') minimal
  surface*, Chem. Phys. Lett. **314** (1999) 543–551.
  [doi:10.1016/S0009-2614(99)01000-3](https://doi.org/10.1016/S0009-2614(99)01000-3)

All three surfaces are Bonnet associates sharing one Weierstrass function

```
R(τ) = 1 / √(τ⁸ − 14τ⁴ + 1)
```

with Bonnet angles `θ_D = 0°`, `θ_G = arccot(K′/K) ≈ 38.0147740°`,
`θ_P = 90°` (`K = K(m=1/4)`, `K′ = K(m=3/4)`, computed by the AGM). The
surface is

```
r(ω) = Re [ e^{iθ} ∫₀^ω (1−τ², i(1+τ²), 2τ) R(τ) dτ ]
```

over the fundamental domain bounded by the two coordinate axes and the
arc `|ω + (1+i)/√2| = √2`, with corners `O = 0`, `P = (√3−1)/√2`,
`R = iP` (P, R are branch points of `R(τ)` — flat points of the surface).

Implementation details:

- The curved-triangle domain is treated as a curved quad `O–P–Q–R`
  (Q = arc midpoint) and mapped from the unit square by a **Coons patch**,
  so the regular `(u,v)` grid becomes the quad mesh of the surface.
- Integrals use composite 16-point **Gauss–Legendre** quadrature marched
  across the grid. Segments near branch points are integrated in the
  local chart `ζ = √(τ − branch)` where the integrand is analytic — full
  accuracy at the singular corners.
- The isometries assembling the patch into one cubic cell were **derived
  numerically from first principles**: the surface was analytically
  continued across the three patch boundary curves; the side-pairing
  isometries were extracted by Procrustes fits (residual ~1e-9); the
  space group was generated by composition; the cubic lattice was read
  off from its pure translations. (The papers' own assembly tables rest
  on internally inconsistent frame descriptions — e.g. the Gyroid's true
  cubic period is 4a in the paper's notation, not 2a.) This reproduces
  the expected crystallography exactly:

  | surface | space group | patches / cell | lattice |
  |---|---|---|---|
  | Gyroid | *Ia3̄d* (230) | 96 | bcc |
  | Schwarz P | *Im3̄m* (229) | 48 | bcc |
  | Schwarz D | *Fd3̄m* (227) | 192 | fcc |

  For P and D the fundamental patch carries an internal 2-fold symmetry
  (it spans two asymmetric units); the pipeline detects the stabilizer
  and quotients it out, so no face is emitted twice. Each operation also
  carries an orientation sign — whether it swaps the two sides of the
  surface — which for P/D is *not* the determinant (their groups contain
  in-surface 2-fold axes and surface-orthogonal mirrors), so it is
  measured during derivation and stored in the table.
- The cubic lattice constants come out as exact elliptic-integral
  expressions (matching the integration to ~1e-9): `L_P = K′`,
  `L_D = 2K`, `L_G = 2KK′ / √(K² + K′²)` in κ=1 integration units.
- Each cell is aligned with its standard level-set convention
  (coordinates in `2π / cell` units), up to the true deviation of those
  nodal approximations: Gyroid `sin x cos y + sin y cos z + sin z cos x`
  (mean |F| ≈ 0.007), P `cos x + cos y + cos z` (≈ 0.045), D
  `sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z`
  (≈ 0.002).

Verified properties of the emitted mesh:

| check | result (all three surfaces) |
|---|---|
| faces | 100 % quads, zero degenerate |
| orientation | globally consistent — on the 3-torus every directed edge appears exactly once |
| minimality | numeric mean curvature ≈ 0 (finite-difference check) |
| periodicity | on the 3-torus every edge borders exactly 2 quads (no gaps, no overlaps) |
| seams | all interior patch boundaries weld exactly (~1e-9) |

### Files

- `weierstrass.py` — pure-numpy core: Weierstrass integration, Coons
  domain, per-surface space-group tables, unit-cell assembly. Importable
  and testable outside Blender.
- `operators.py` — `tpms.generate`: builds the mesh object + Array stack.
- `properties.py`, `ui.py` — settings and N-panel.

Schoen I-WP and F-RD have no published exact parametrization from this
series (only nodal approximations); adding them would require deriving
Weierstrass data with cube-root branch points — possible future work on
the same footing.

---

## 中文

三种曲面都由精确的 Weierstrass 表示直接求值：每个顶点都严格落在真正的
极小曲面上（平均曲率恒为零），四边形沿曲面天然参数线分布，网格精确周期，
平铺的晶胞无缝焊合。常见的等值面近似 `sin x cos y + sin y cos z +
sin z cos x = 0` 在这里只作为验证参照出现。

### 用法

1. 安装插件（Edit → Preferences → Add-ons → Install from Disk）。
2. 3D 视口 N 面板 → **TPMS** 标签页。
3. 选择*类型*（Gyroid / Schwarz P / Schwarz D），设置*晶胞尺寸*、
   *X/Y/Z 晶胞数*、*分辨率* → **Generate TPMS**。

结果是一个晶胞网格加三个 Array 修改器。晶胞数量可以在修改器上实时修改；
需要单一网格时应用修改器即可。晶胞边界的周期性精确到晶胞尺寸的 ~1e-9，
所以平铺后的接缝无缝且水密。

- **分辨率** —— 每个基本片边上的四边形数。一个晶胞 = `片数 × res²` 个四
  边形，Gyroid / P / D 分别为 96 / 48 / 192 片（res 8 → 6144 / 3072 /
  12288）。任何分辨率下顶点都是精确的；提高 res 只是让轮廓更平滑。
- **平滑着色** —— 使用精确解析法向（Weierstrass 数据的 Gauss 映射）。

关于块体边缘的说明：晶胞是由**完整**基本片构成的一个精确平移单元。
Schwarz D 的基本片恰好终止在晶胞面上，Gyroid 的外伸
可以忽略，所以它们的块体看起来是方正的。Schwarz P 的基本片会斜跨晶胞面
（P 曲面内嵌有 `(t, t+1/2, 1/4)·L` 这类直线，穿透晶胞面），因此有限大小
的 P 块体表皮是参差的——尽管平铺本身完全无缝。需要平切块体时用 Boolean
修整即可。

### 数学原理

基于 Klinowski 等人的精确计算系列论文：

- **Gyroid** — P.J.F. Gandy, J. Klinowski, *Exact computation of the
  triply periodic G ('Gyroid') minimal surface*, Chem. Phys. Lett.
  **321** (2000) 363–371.
  [doi:10.1016/S0009-2614(00)00373-0](https://doi.org/10.1016/S0009-2614(00)00373-0)
  （[PDF](https://mathcurve.com/surfaces.gb/Gyroide/sdarticle%20gyroid.pdf)）
- **Schwarz P** — P.J.F. Gandy, J. Klinowski, *Exact computation of the
  triply periodic Schwarz P minimal surface*, Chem. Phys. Lett. **322**
  (2000) 579–586.
  [doi:10.1016/S0009-2614(00)00453-X](https://doi.org/10.1016/S0009-2614(00)00453-X)
  （[PDF](https://www.mathcurve.com/surfaces.gb/schwarz/sdarticle.pdf)）
- **Schwarz D** — P.J.F. Gandy, D. Cvijović, A.L. Mackay, J. Klinowski,
  *Exact computation of the triply periodic D ('diamond') minimal
  surface*, Chem. Phys. Lett. **314** (1999) 543–551.
  [doi:10.1016/S0009-2614(99)01000-3](https://doi.org/10.1016/S0009-2614(99)01000-3)

三种曲面是 Bonnet 伴随家族，共享同一个 Weierstrass 函数

```
R(τ) = 1 / √(τ⁸ − 14τ⁴ + 1)
```

Bonnet 角分别为 `θ_D = 0°`、`θ_G = arccot(K′/K) ≈ 38.0147740°`、
`θ_P = 90°`（`K = K(m=1/4)`、`K′ = K(m=3/4)`，用 AGM 计算）。曲面为

```
r(ω) = Re [ e^{iθ} ∫₀^ω (1−τ², i(1+τ²), 2τ) R(τ) dτ ]
```

基本域由两条坐标轴和圆弧 `|ω + (1+i)/√2| = √2` 围成，角点 `O = 0`、
`P = (√3−1)/√2`、`R = iP`（P、R 是 `R(τ)` 的分支点——曲面的平点）。

实现细节：

- 曲边三角域被视为曲边四边形 `O–P–Q–R`（Q 为弧中点），用 **Coons patch**
  从单位正方形映射，正方形的规则 `(u,v)` 网格直接成为曲面的四边形网格。
- 积分采用逐段推进的 16 点 **Gauss–Legendre** 复合求积。靠近分支点的线段
  换到局部坐标 `ζ = √(τ − 分支点)` 中积分（被积函数在该坐标下解析），奇
  异角点处不损失精度。
- 把基本片拼成一个立方晶胞的等距变换是**从第一性原理数值推导的**：把
  曲面沿基本片的三条边界曲线解析延拓，用 Procrustes 拟合提取边配对等距
  （残差 ~1e-9），复合生成空间群，再从纯平移读出立方晶格。（论文自己的
  拼装表格建立在内部矛盾的坐标系描述上——例如 Gyroid 的真实立方周期是论
  文记号下的 4a，而非 2a。）结果与教科书晶体学完全吻合：

  | 曲面 | 空间群 | 片数/晶胞 | 晶格 |
  |---|---|---|---|
  | Gyroid | *Ia3̄d* (230) | 96 | 体心 |
  | Schwarz P | *Im3̄m* (229) | 48 | 体心 |
  | Schwarz D | *Fd3̄m* (227) | 192 | 面心 |

  P 和 D 的基本片带有内部 2 重对称（跨越两个不对称单元）；管线会检测这
  个稳定子群并做商，保证没有面被重复输出。每个操作还带有一个定向符号
  ——它是否交换曲面两侧——对 P/D 这**不等于**行列式（它们的对称群包含躺
  在曲面内的 2 重轴和垂直于曲面的镜面），所以在推导阶段逐一测定并存入
  表中。
- 立方晶格常数有精确的椭圆积分表达式（与数值积分吻合到 ~1e-9）：
  `L_P = K′`、`L_D = 2K`、`L_G = 2KK′ / √(K² + K′²)`（κ=1 积分单位）。
- 每个晶胞都与其标准等值面约定对齐（坐标以 `2π / 晶胞` 为单位），偏差即
  这些 nodal 近似的真实误差：Gyroid `sin x cos y + sin y cos z +
  sin z cos x`（mean |F| ≈ 0.007）、P `cos x + cos y + cos z`（≈ 0.045）、
  D `sin x sin y sin z + sin x cos y cos z + cos x sin y cos z +
  cos x cos y sin z`（≈ 0.002）。

输出网格的验证结果：

| 检查项 | 结果（三种曲面均通过） |
|---|---|
| 面 | 100% 四边形，零退化 |
| 定向 | 全局一致——3-环面上每条有向边恰好出现一次 |
| 极小性 | 数值平均曲率 ≈ 0（有限差分检验） |
| 周期性 | 3-环面上每条边恰好邻接 2 个面（无缺口、无重叠） |
| 接缝 | 所有内部基本片边界精确焊合（~1e-9） |

### 文件

- `weierstrass.py` —— 纯 numpy 核心：Weierstrass 积分、Coons 参数域、
  各曲面的空间群表、晶胞组装。可在 Blender 外独立导入和测试。
- `operators.py` —— `tpms.generate`：构建网格对象和 Array 修改器栈。
- `properties.py`、`ui.py` —— 参数设置和 N 面板。

Schoen I-WP 和 F-RD 在该系列中没有发表过精确参数化（只有 nodal 近似）；
要加入它们需要推导带立方根分支点的 Weierstrass 数据——可作为同一思路下
的后续工作。

## License

MIT
