# PureQuad TPMS

[English README](README.md)

面向 Blender 的原生四边拓扑 TPMS：解析生成、网格轻量、方便编辑。

![Gyroid、Schwarz P 与 Schwarz D 全四边网格](docs/tpms_showcase.png)

PureQuad TPMS 根据精确的 Enneper–Weierstrass 参数化直接生成 Gyroid、
Schwarz P 和 Schwarz D。它不从体素场提取等值面再做清理，而是从一开始
就按照曲面自身的结构构造全四边网格。

## PureQuad 的优势

### 简单、干净的拓扑

所有面都是四边形，边流沿曲面自身的参数线排列。精确基本三角片两两配对，
形成真正的四边宏观片：Gyroid 48 张、Schwarz P 24 张、Schwarz D 96 张。
配对三角片的公共圆弧位于宏观片内部，不会成为 Blender 中的拓扑边。

生成结果轻量而且结构确定，没有体素噪声、任意三角剖分，也不需要重拓扑。
可以直接用于循环边选择、Subdivision Surface、Solidify、变形和其他常规
Blender 网格操作。

### 精确曲面求值

每个顶点都由极小曲面的精确 Enneper–Weierstrass 表示求得。降低网格密度
只会简化拓扑，不会把曲面替换成粗糙的体素近似。

### 快速生成

PureQuad TPMS 采样二维解析曲面片，不需要建立和提取三维体素场。在常用分
辨率下可以即时生成。插件为纯 Python，只使用 Blender 自带的 NumPy，没有
外部依赖。

### 无缝周期平铺

生成器创建一个立方晶胞，并自动附加 X、Y、Z 三个方向的 Array 修改器。
周期边界顶点精确对应，相邻晶胞之间没有可见接缝，可以合并成连续点阵。

## 安装

1. 从[最新 Release](https://github.com/Zhangyanbo/purequad-tpms/releases/latest)下载 ZIP。
2. 在 Blender 中打开 **Edit → Preferences → Add-ons → 菜单 → Install from Disk**。
3. 选择下载的 ZIP 并启用扩展。

需要 Blender 4.2 或更高版本。

## 使用

在 3D 视口打开 N 面板，选择 **TPMS** 标签页。设置曲面类型、晶胞尺寸、平
铺数量和 **Quad Subdivisions**，然后点击 **Generate TPMS**。

Quad Subdivisions 表示每张宏观片每边的四边面数量。细分值为 `n` 时，每个
晶胞的面数为：

- Gyroid：`48 × n²`
- Schwarz P：`24 × n²`
- Schwarz D：`96 × n²`

提高细分值会使轮廓更加平滑，但宏观拓扑保持不变。Gyroid 和 Schwarz D 可
以形成边界方正的有限块体。Schwarz P 的曲面片会跨过晶胞边界，因此有限块
体的外表面天然错落，但周期平铺本身仍然无缝；需要平整外表面时可以使用
Boolean 切割。

### 等值面偏移

**Iso Level** 把曲面移动到标准三角函数 TPMS 场的等值面 `F(x, y, z) = t`
上，从而调节两侧通道的体积比例——即体素类 TPMS 生成器中常见的
isolevel 参数。默认 `t = 0` 时网格仍然位于精确极小曲面上；其他取值沿场
的梯度搬运每个顶点，因此全四边网格、宏观片布局和无缝平铺完全不变。取值
范围到颈部收缩为止：Gyroid 约 `|t| < 1.3`，Schwarz P 和 D 约
`|t| < 0.99`。

## 数学背景

曲面依据 Gandy、Cvijović、Mackay 与 Klinowski 的精确计算工作。每张真正
的宏观四边片都由显式解析方形映射描述，并通过圆反演连接成对的基本域。晶
胞等距变换通过数值方法推导，并依据相应空间群进行验证。

完整公式、构造方法和验证细节见[数学文档](docs/mathematics.md)。

## 许可证

MIT
