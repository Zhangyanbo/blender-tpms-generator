# TPMS Generator

> Blender 4.2+ extension that fills any target mesh with a solid Triply
> Periodic Minimal Surface (Gyroid, Schwarz P, Schwarz D, Schoen IWP,
> Fischer-Koch S).

A Blender 4.2+ 扩展插件，可在用户选定的目标网格内部，生成五种类型 TPMS
（默认 Gyroid）的实体网格，参数可调。适用于轻量化结构、3D 打印晶格、
仿生散热、骨支架等场景。

![preview placeholder](docs/preview.png)

---

## Features

- **5 TPMS families** — Gyroid (default), Schwarz Primitive, Schwarz Diamond,
  Schoen I-WP, Fischer-Koch S — all defined by their canonical level-set
  equations.
- **Two solid modes**:
  - *Shell* — uniform wall thickness around `f = iso` (good for 3D printing,
    heat exchangers).
  - *Volume* — one phase (`f < iso`) as a solid (good for closed-cell
    lattices).
- **Field-domain clipping to the target** — combines the TPMS scalar field
  with the target's signed distance field (`g = max(g_tpms, sd)`) instead
  of a Boolean modifier. Works on non-closed targets like default Suzanne
  (open neck, eye sockets, floating earring rings).
- **5-ray majority in/out test** — robust against open caps and small holes.
- **Two-pass SDF computation** — coarse sub-sampled grid + narrow-band
  refinement at the iso-boundary; tractable for million-voxel grids.
- **Analytic vertex projection** — after Naive Surface Nets extracts the
  topology, each interior vertex is Newton-stepped onto the true analytic
  iso-surface, giving clean geometry without the wobble of plain Surface
  Nets.
- **Connected-component filter** — drops boundary specks regardless of
  where they came from.
- **Pure Python + numpy** — no external dependencies, no Marching Cubes
  lookup tables, no OpenVDB.

---

## Installation

### From a local folder (recommended for now)

1. Clone or download this repository to a known location.
2. In Blender: **Edit → Preferences → Extensions → ▼ (top-right) →
   Install from Disk…**
3. Select the project folder (the one containing `blender_manifest.toml`)
   or zip it first and pick the `.zip`.
4. Enable **TPMS Generator** in the extensions list.

### Updating after editing the source

Blender will **not** automatically re-import a running add-on. After
editing any `.py` file:

- Disable the add-on (click ✓ to ✗), then re-enable it, **or**
- Restart Blender.

---

## Usage

1. Select (or create) a closed-ish mesh to use as the *container* — for
   example a Suzanne head, an imported STL, or a primitive cube.
2. Open the 3D View → **N panel → TPMS tab**.
3. Set **Target Mesh** to your container.
4. Pick a **Type** (default Gyroid).
5. Tune **Cell Size**, **Iso Value**, **Wall Thickness**.
6. Click **Generate TPMS**.

The new mesh `TPMS_<type>` is added to the active collection.

---

## Parameters

### Lattice

| Parameter      | What it controls |
| -------------- | ---------------- |
| **Cell Size**  | World-space length of one TPMS period (one `2π` cycle). |
| **Iso Value**  | Level-set constant `c`. `0` ≈ equal phase volumes; ±values shift the volume fraction. Range roughly ±1.2. |
| **Origin**     | Phase-shift the lattice. Useful for tiling alignment. |
| **Rotation**   | Rotate the lattice independently of the target. |

### Solid

| Parameter         | What it controls |
| ----------------- | ---------------- |
| **Mode**          | *Shell* — `\|f − iso\| − thickness` (closed double sheet). *Volume* — one phase (`f < iso`) of the TPMS. |
| **Wall Thickness**| For Shell mode. In *field units* (Schwarz P has range ±3, Gyroid ±1.5, etc.). Try 0.2–0.5 for thin lattice walls. |

### Sampling

| Parameter             | What it controls |
| --------------------- | ---------------- |
| **Resolution / Cell** | Voxels per TPMS period along each axis. 24–48 is a good range. |
| **Surface Snap Iters**| Newton iterations that project each Surface-Nets vertex onto the analytic iso-surface. 2–3 is plenty. `0` disables (raw Surface Nets, blurry). |
| **BBox Padding**      | Extra space around the target's bounding box before sampling. |

### Output

| Parameter              | What it controls |
| ---------------------- | ---------------- |
| **Clip to Target**     | If on, the TPMS is restricted to the target's interior via the SDF combination. If off, the raw TPMS is returned with a bbox cap. |
| **Robust In/Out Test** | Use 5-ray majority-vote parity for inside/outside. Slower but mandatory for concave or non-closed targets. |
| **SDF Subsample**      | Compute the target's SDF on every Nth voxel, then refine the narrow boundary band. `4` is a good default; lower for crisper boundary. |
| **Boundary Inset**     | Shrink the target by N voxels before clipping. Use 0.5–1.5 to hide any residual fuzz. |
| **Min Component Faces**| Drop isolated mesh fragments smaller than this many quads. 200 is the default; the legitimate TPMS body is always 10⁴–10⁵ quads. |
| **Smooth Shading**     | Apply smooth shading to the generated mesh. |

---

## How it works

```
target mesh
    │
    │  BVH (mathutils.bvhtree)
    ▼
signed distance field  sd(x,y,z)         analytic TPMS field  f(x,y,z)
   (5-ray majority in/out,                 (sin/cos formulas,
    coarse + boundary refinement)           per-type scalar field)
    │                                       │
    │            shell:  g_tpms = |f − iso| − thickness
    │            volume: g_tpms = f − iso
    │                                       │
    └──────────►  g = max(g_tpms, sd)  ◄────┘
                       │
                       │  Naive Surface Nets (custom)
                       ▼
                  quad mesh of g = 0
                       │
                       │  Newton iterations against ∇f
                       ▼
                vertices snapped onto the analytic iso-surface
                       │
                       │  Union-find connected components
                       ▼
                final manifold quad mesh
```

The TPMS surfaces are implicit (level-set), not parametric — there's no
closed-form `(u,v) → ℝ³` mapping. So the pipeline is fundamentally
"sample a scalar field on a regular grid, extract an iso-surface". Mesh
quality therefore scales with the sampling resolution and the cleanliness
of the iso-surface extractor.

The choice of Naive Surface Nets + analytic vertex projection (over
Marching Cubes) gives:

- Similar surface accuracy to MC at the saddles that dominate TPMS
  topology.
- A ~10× shorter implementation (no 256×16 triangle table).
- Quad output (better than MC's triangles for downstream subdivision /
  remeshing).

Field-domain clipping (instead of a `Boolean INTERSECT` modifier) is what
makes the plugin survive on default Suzanne and other "almost-closed"
meshes that EXACT Boolean falls apart on.

---

## TPMS equations

All five surfaces are written here in lattice coordinates (one period =
`2π`):

```
Gyroid          : sin x · cos y + sin y · cos z + sin z · cos x
Schwarz P       : cos x + cos y + cos z
Schwarz D       : sin x sin y sin z + sin x cos y cos z
                + cos x sin y cos z + cos x cos y sin z
Schoen I-WP     : 2 (cos x cos y + cos y cos z + cos z cos x)
                − (cos 2x + cos 2y + cos 2z)
Fischer-Koch S  : cos 2x sin y cos z + cos x cos 2y sin z
                + sin x cos y cos 2z
```

---

## Limitations / known issues

- The BVH SDF pass is a Python loop — large meshes × high resolution can
  take 10–30 s. If it feels slow, drop **Resolution / Cell** to 28–32
  while you're tuning, then crank it up for the final render.
- Open targets with very thin walls (single-sided sheets) can still
  confuse the in/out test in places; the component filter is the safety
  net.
- The mesh is closed and manifold for downstream Boolean / Solidify /
  Subdivide use, but is **not** guaranteed to have well-formed normals
  across all saddles. Apply *Shade Smooth* and (optionally) *Mesh →
  Normals → Recalculate Outside* if your renderer cares.

---

## License

[MIT](LICENSE).

> ⚠️ Note: the Blender Foundation's legal FAQ argues that add-ons calling
> the `bpy` Python API are derivative works of Blender and should
> therefore be GPL-2.0-or-later. This is a long-standing grey area for
> Python add-ons and is rarely enforced in practice. If you redistribute
> this add-on alongside Blender you may want to switch the license to
> GPL-2.0-or-later to stay on the conservative side.

---

## Acknowledgements

- TPMS equations: Schoen (NASA TN D-5541, 1970), Schwarz, Fischer-Koch.
- Surface Nets: S. F. Gibson, *Constrained elastic surface nets*, MERL TR
  98-19 (1998).
