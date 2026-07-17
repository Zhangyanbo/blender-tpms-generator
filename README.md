# TPMS Generator

**Exact** triply periodic minimal surfaces for Blender 4.2+ — Gyroid,
Schwarz P and Schwarz D as clean, all-quad meshes from their analytic
Enneper–Weierstrass parametrizations.

Blender 4.2+ 扩展插件 —— 由解析 Enneper–Weierstrass 参数化生成的**精确**
三周期极小曲面（Gyroid、Schwarz P、Schwarz D），输出干净的全四边形网格。

| Gyroid (2×2×2) | Schwarz P (2×2×2) | Schwarz D (2×2×2) |
|:---:|:---:|:---:|
| ![Gyroid](docs/gyroid_2x2x2.png) | ![Schwarz P](docs/schwarz_p_2x2x2.png) | ![Schwarz D](docs/schwarz_d_2x2x2.png) |

| Gyroid cell | Schwarz P cell | Schwarz D cell |
|:---:|:---:|:---:|
| ![Gyroid quad flow](docs/gyroid_quadflow.png) | ![Schwarz P cell](docs/schwarz_p_cell.png) | ![Schwarz D cell](docs/schwarz_d_cell.png) |

[English](#english) | [中文](#中文)

---

## English

TPMS meshes usually come from running marching cubes over a level-set
approximation, which produces dense triangle meshes that only sit near
the true surface. This add-on takes the analytic route instead: it
integrates the Enneper–Weierstrass representation of each surface
directly, so the mesh you get is the mathematical object itself.

Concretely, that means:

- every vertex lies on the true minimal surface — mean curvature is
  exactly zero at any resolution;
- faces are 100 % quads following the surface's natural parameter
  lines, ready for subdivision and editing with no retopology;
- the unit cell is exactly periodic, so tiling it with ordinary Array
  modifiers gives seamless, watertight joints;
- shading uses exact analytic normals from the Gauss map, independent
  of mesh density;
- meshes stay light — a faithful Gyroid cell is ~6k quads and builds in
  milliseconds.

The mathematics behind all of this — the Weierstrass representation,
the numerically derived space-group assembly, and the verification
suite — is written up in [docs/mathematics.md](docs/mathematics.md),
based on the exact-computation papers by Gandy, Cvijović, Mackay &
Klinowski.

### Setting it up

1. Download this repository (Code → Download ZIP, or `git clone`).
2. In Blender: Edit → Preferences → Add-ons → ⌄ → *Install from Disk*,
   and select the ZIP / folder.
3. Enable **TPMS Generator**.

Blender 4.2+ is required. The add-on is pure Python + numpy (bundled
with Blender), with no external dependencies.

### Generating a lattice

Open the N-panel in the 3D viewport and switch to the **TPMS** tab.
Pick a surface type, set the cell size, the number of cells along each
axis, and the resolution, then press **Generate TPMS**. The operator
creates one unit-cell mesh with three Array modifiers on it — cell
counts stay editable on the modifiers, and applying them collapses the
lattice into a single mesh.

A few things worth knowing:

- *Resolution* is the number of quads per fundamental-patch edge; a
  cell consists of 96 / 48 / 192 patches for Gyroid / P / D. Since
  vertices are exact at any resolution, raising it only smooths
  silhouettes.
- The cell is a translational unit built from whole fundamental
  patches. Gyroid and Schwarz D blocks come out box-clean; Schwarz P
  patches straddle the cell faces, so a finite P block has a ragged
  skin even though its tiling is exactly seamless. Trim with a Boolean
  when you need a flat-cut P block.

---

## 中文

常见的 TPMS 网格来自对等值面近似做 marching cubes，得到的是只贴近真实
曲面的稠密三角网格。本插件走解析路线：直接积分每种曲面的
Enneper–Weierstrass 表示，得到的网格就是数学对象本身。

具体来说：

- 每个顶点都严格落在真正的极小曲面上——任何分辨率下平均曲率都恒为零；
- 面 100% 是四边形，沿曲面天然参数线分布，可直接细分和编辑，无需重拓
  扑；
- 晶胞精确周期，用普通 Array 修改器平铺即可得到无缝、水密的接合；
- 着色使用来自 Gauss 映射的精确解析法向，与网格密度无关；
- 网格轻量——一个忠实的 Gyroid 晶胞约 6k 四边形，毫秒级生成。

这一切背后的数学——Weierstrass 表示、数值推导的空间群拼装、验证套件
——整理在 [docs/mathematics.md](docs/mathematics.md) 中，基于 Gandy、
Cvijović、Mackay 与 Klinowski 的精确计算系列论文。

### 装好它

1. 下载本仓库（Code → Download ZIP，或 `git clone`）。
2. 在 Blender 中：Edit → Preferences → Add-ons → ⌄ → *Install from
   Disk*，选择 ZIP / 文件夹。
3. 启用 **TPMS Generator**。

需要 Blender 4.2+。插件为纯 Python + numpy（Blender 自带），无外部依
赖。

### 生成点阵

在 3D 视口打开 N 面板，切换到 **TPMS** 标签页。选择曲面类型，设置晶胞
尺寸、各方向晶胞数和分辨率，点击 **Generate TPMS**。操作符会创建一个
挂着三个 Array 修改器的晶胞网格——晶胞数量可以在修改器上继续改，应用修
改器即可把点阵合并成单一网格。

几点值得了解：

- *分辨率*是每个基本片边上的四边形数；一个晶胞由 96 / 48 / 192 个基本
  片组成（Gyroid / P / D）。顶点在任何分辨率下都是精确的，提高分辨率只
  是让轮廓更平滑。
- 晶胞是由完整基本片构成的平移单元。Gyroid 和 Schwarz D 的块体是方正
  的；Schwarz P 的基本片会斜跨晶胞面，因此有限大小的 P 块体表皮参差
  （平铺本身仍完全无缝）。需要平切的 P 块体时用 Boolean 修整。

## License

MIT
