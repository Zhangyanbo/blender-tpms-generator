# Repository Working Agreement

These rules apply to every agent changing the PureQuad TPMS repository.

## Definition of done

A plugin code change is not finished when the source files or unit tests are
finished. It is finished only after all of the following have happened:

1. Run the smallest relevant source-level test suite.
2. Build the local Blender extension package.
3. Install and enable that package in the user's installed Blender.
4. Load the plugin from Blender's actual user-extension directory and perform
   one minimal generation check there.

Installation is the default final step for every plugin change. Do not pause
to ask whether the user wants installation. OS or sandbox approval may still
be requested when the tool requires it. Commit, tag, push, and GitHub Release
are different: perform them only after explicit user approval.

## Find and use the real Blender installation

- The Blender executable does not need to be inside this repository. On this
  machine, check `/Applications/Blender.app/Contents/MacOS/Blender` and other
  normal application locations before claiming that Blender is unavailable.
- Detect the installed Blender version and use its actual API identifiers.
  Do not assume that an enum or property name from another Blender version is
  valid. Prefer one small introspection check over repeated failed launches.
- Prefer Blender's official extension build and install commands. Verify the
  installed manifest version and, where useful, compare installed source files
  with the workspace before testing the installed copy.
- Consolidate validation into one background Blender process whenever possible.

## Geometry and numerical method

- Start from the mathematical structure of the surface. For a genuine macro
  quad, first seek one direct analytic expression over the whole quadrilateral.
- Do not route a quad through a triangle mesh, harmonic reparameterization,
  polynomial fit, lookup table, or inverse solve when an exact analytic map is
  available.
- Do not expose implementation artifacts such as solver resolution or
  quadrature order in the UI unless users genuinely need to control them.
- Gyroid, Schwarz P, and Schwarz D are Bonnet associates. Reuse their common
  macro-domain construction instead of maintaining three unrelated pipelines.
- Treat interaction latency as a product requirement. A result under one
  second can still be too slow if an earlier button press felt instantaneous.
  Measure cold and warm generation time, and keep the click-to-result path
  direct.

## Testing and experimentation

- Use focused numerical assertions: topology, all-quad faces, periodic seams,
  positive parameter-domain Jacobian, analytic branch continuity, and Blender
  object creation.
- Do not start Matplotlib, font-cache generation, benchmark plotting, or a
  separate visualization pipeline unless the user explicitly asks for a plot.
  Timings should normally be printed as numbers.
- Do not create an experiment branch, temporary worktree, or large benchmark
  harness when a small direct check answers the question.
- Remove temporary scripts, worktrees, caches, and intermediate outputs created
  during the task.

## README imagery

- Keep the README focused on one strong composite image rather than a gallery
  of redundant views.
- When showcasing this plugin, render with Blender and show the actual generated
  geometry. The visible wire lines must correspond to real mesh edges.
- The primary visual message is the clean, regular all-quad topology. Avoid
  labels and decorative text inside the image when the geometry is sufficient.
- Visually inspect the final render at full resolution. Reject clipped objects,
  missing wires, weak contrast, excessive empty space, or a composition that
  makes the mesh secondary.

## Workflow and communication

- Follow the user's requested order of operations. If asked to release the
  current code before starting a documentation change, keep those histories
  separate.
- Do not stop while a safe, required workflow step remains. In particular, do
  not confuse the need for commit/release approval with the standing requirement
  to install every plugin update.
- Report the exact state: source tested, locally packaged, installed, verified
  from the installed copy, committed, or released. These are distinct states.
- Prefer one well-reasoned direct implementation over a chain of approximations
  whose complexity is paid back as runtime cost and maintenance burden.
