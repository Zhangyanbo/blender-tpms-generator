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

### Why

Most TPMS tools mesh a level-set approximation with marching cubes and
hand you a large triangle soup. This add-on evaluates the surfaces from
their exact Weierstrass representations instead, which buys you:

- **Every vertex on the true minimal surface** — mean curvature is
  exactly zero, at any resolution.
- **All-quad meshes along the surface's natural parameter lines** —
  clean topology for subdivision, editing, retopo-free workflows.
- **Exact periodicity** — one unit cell tiled by ordinary Array
  modifiers welds seamlessly and watertight across cell boundaries.
- **Exact analytic normals** — smooth shading from the Gauss map of the
  surface, independent of mesh density.
- **Light meshes** — a faithful Gyroid cell is ~6k quads at the default
  resolution and builds in milliseconds.

The mathematics behind the add-on (Weierstrass representation, numerical
derivation of the space-group assembly, verification) is documented in
[docs/mathematics.md](docs/mathematics.md), based on the exact-computation
papers by Gandy, Cvijović, Mackay & Klinowski.

### Install

1. Download this repository (Code → Download ZIP, or `git clone`).
2. Blender → Edit → Preferences → Add-ons → ⌄ → *Install from Disk*,
   select the folder / ZIP.
3. Enable **TPMS Generator**.

Requires Blender 4.2+. Pure Python + numpy (bundled with Blender), no
external dependencies.

### Usage

1. Open the N-panel in the 3D viewport → **TPMS** tab.
2. Choose *Type* (Gyroid / Schwarz P / Schwarz D), set *Cell Size*,
   *Cells X/Y/Z* and *Resolution*.
3. Press **Generate TPMS**.

You get one unit-cell mesh plus three Array modifiers. Cell counts are
editable live on the modifiers; apply them when you want a single mesh.

- **Resolution** is the quad count per fundamental-patch edge; one cell
  has 96 / 48 / 192 patches for Gyroid / P / D. Vertices are exact at
  any resolution — raise it only for smoother silhouettes.
- The cell is a translational unit built from whole fundamental patches.
  Gyroid and Schwarz D blocks look box-clean; Schwarz P patches straddle
  the cell faces, so a finite P block has a ragged skin even though the
  tiling itself is exactly seamless — trim with a Boolean for a flat-cut
  block.

---

## 中文

### 为什么用它

多数 TPMS 工具对等值面近似做 marching cubes，得到巨大的三角面汤。本插件
改为直接从精确的 Weierstrass 表示求值曲面，好处是：

- **每个顶点都严格落在真正的极小曲面上** —— 任何分辨率下平均曲率都恒为
  零。
- **沿曲面天然参数线的全四边形网格** —— 拓扑干净，适合细分、编辑，无需
  重拓扑。
- **精确周期性** —— 一个晶胞用普通 Array 修改器平铺，跨晶胞边界无缝且水
  密焊合。
- **精确解析法向** —— 平滑着色来自曲面的 Gauss 映射，与网格密度无关。
- **网格轻量** —— 默认分辨率下一个忠实的 Gyroid 晶胞约 6k 四边形，毫秒
  级生成。

插件背后的数学（Weierstrass 表示、空间群拼装的数值推导、验证）见
[docs/mathematics.md](docs/mathematics.md)，基于 Gandy、Cvijović、
Mackay 与 Klinowski 的精确计算系列论文。

### 安装

1. 下载本仓库（Code → Download ZIP，或 `git clone`）。
2. Blender → Edit → Preferences → Add-ons → ⌄ → *Install from Disk*，
   选择文件夹 / ZIP。
3. 启用 **TPMS Generator**。

需要 Blender 4.2+。纯 Python + numpy（Blender 自带），无外部依赖。

### 用法

1. 3D 视口打开 N 面板 → **TPMS** 标签页。
2. 选择*类型*（Gyroid / Schwarz P / Schwarz D），设置*晶胞尺寸*、
   *X/Y/Z 晶胞数*和*分辨率*。
3. 点击 **Generate TPMS**。

结果是一个晶胞网格加三个 Array 修改器。晶胞数量可以在修改器上实时修改；
需要单一网格时应用修改器即可。

- **分辨率**是每个基本片边上的四边形数；一个晶胞的片数 Gyroid / P / D
  分别为 96 / 48 / 192。任何分辨率下顶点都是精确的——提高分辨率只是让轮
  廓更平滑。
- 晶胞是由完整基本片构成的平移单元。Gyroid 和 Schwarz D 的块体是方正
  的；Schwarz P 的基本片会斜跨晶胞面，因此有限大小的 P 块体表皮参差（平
  铺本身仍完全无缝）——需要平切块体时用 Boolean 修整。

## License

MIT
