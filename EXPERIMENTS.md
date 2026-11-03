# Experiment index

This index maps the internal stage IDs to their experiment topics and code snapshots. The IDs are retained for tracing records; the descriptions below spell out what each stage examines.

| Stage ID | Experiment | Code and records |
|---|---|---|
| S0 / T0 | Neo-Hookean material derivatives and element assembly; implicit-Euler updates for a linear oscillator | [`material_validation.cpp`](simulator/cpp/examples/material_validation.cpp), [`dynamic_reference_check.cpp`](simulator/cpp/examples/dynamic_reference_check.cpp); related records under [`simulator/docs/`](simulator/docs/) |
| S1 | Uniform rod stretching and diagnostics of the Dirichlet boundary condition (DBC) penalty stiffness | [`tension_check.cpp`](simulator/cpp/examples/tension_check.cpp); [diagnostic workspace](experiments/workspaces/fiction-s1-diagnostic-20261031-01/) and [two-group workspace](experiments/workspaces/fiction-s1-two-group-20261031-01/) |
| S2 | Static cantilever bending and spatial mesh refinement | [`cantilever_check.cpp`](simulator/cpp/examples/cantilever_check.cpp); [S2 extension workspace](experiments/workspaces/fiction-s2-extension-20261031-01/) |
| S3 | Dynamic cantilever bending and time-step refinement | [Original S3 workspace](experiments/workspaces/fiction-s3-20261031-01/), [S3 time-step and mesh workspace](experiments/workspaces/fiction-s3-dtmesh-20261102-01/), [S3 algebraic-accuracy workspace](experiments/workspaces/fiction-s3-algebra-main-20261102-01/) |
| S3t | Prescribed twisting motion and visualization drafts | [S3/S3t visual workspace](experiments/workspaces/fiction-s3-s3t-visual-main-20261102-01/) |
| S4 | Deformation of a non-rod torus geometry | [S4 workspace](experiments/workspaces/fiction-s4-20261031-01/) |
| S5 | Armadillo mesh scale pilot and larger-mesh run | [S5 pilot workspace](experiments/workspaces/fiction-s5-20261031-01/), [54K main-run workspace](experiments/workspaces/fiction-s5-54k-main-20261102-01/) |
| S6a | Contact-free one-body, four-body, and 16-body scenes | [Original S6a workspace](experiments/workspaces/fiction-s6a-20261031-01/), [single-body](experiments/workspaces/fiction-s6a-single-main-20261102-01/), [four-body](experiments/workspaces/fiction-s6a-four-main-20261102-01/), and [16-body](experiments/workspaces/fiction-s6a-sixteen-main-20261102-01/) workspaces |
| S6b | Contact investigation in isolated workspaces | [C0 workspace](experiments/workspaces/fiction-s6b-c0-20261102-01/), [C1 workspace](experiments/workspaces/fiction-s6b-c1-20261102-01/); the external IPC Toolkit is not included |
| Stopping criteria | Isolated diagnostics of dynamic physical force residuals and solver stopping rules | [Initial diagnostics](experiments/workspaces/fiction-stopcriterion-20261031-01/), [corrected diagnostics](experiments/workspaces/fiction-stopcriterion-corrected-20261031-01/), and [combined-criterion comparisons](experiments/workspaces/fiction-s2s6a-combined-20261031-01/) |

The [`experiments/workspaces/`](experiments/workspaces/) snapshots include selected program files, headers, build configuration, scripts, lightweight text files, and small test meshes. They do not include the workspaces' build directories, binaries, third-party checkouts, or original results. Historical scripts may rely on machine-specific paths and large meshes that are not included here. The original full data, source archives, and experiment reports remain under `simulator/output/validation/` in the original local fiction workspace; this repository does not replace them.

## Large meshes not included in this repository

| File | SHA-256 | Historical input location on the original machine |
|---|---|---|
| Armadillo13K.msh | `03beaf3d4c9692a10664a3a90f1baa8c10d43d7ba47126437d420b68a36e6bce` | `fiction/assets/ipc/tetMeshes/Armadillo13K.msh` |
| Armadillo28K.msh | `73395043488c6c1a99117ce750e4a7ac1f2f3d8849ccd8f18be277a99c0e359b` | `fiction/assets/ipc/tetMeshes/Armadillo28K.msh` |
| Armadillo54K.msh | `e1faa13af8543e01bda3620ae83fe9a1e5ddff1280861d5896b850a84493d520` | `Archive/IPC/input/tetMeshes/Armadillo54K.msh` |

The corresponding scenes reuse the mesh coordinates and prescribed loading directions from the original study. When running them, check the applicable specification under `simulator/docs/` and the experiment script.
