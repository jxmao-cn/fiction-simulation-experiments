# fiction-simulation-experiments

Source code and experiment workspaces for fiction, a 3D solid simulator. The repository covers material and time-integration checks, uniform stretch, cantilever bending, twisting, non-rod geometry, scale experiments, and contact-free multi-body scenes. It is currently a private research snapshot, organized so readers can inspect the implementations and configurations used in the documented experiments.

## Start here

| Read or browse | Contents |
|---|---|
| [Experiment article (English)](blog/fiction-en.md) · [中文](blog/fiction.md) | Narrative overview with figures; assets are shared under [`blog/assets/`](blog/assets/) |
| [Blog file guide](blog/README.md) | Article, asset, editable-figure, and source organization |
| [Experiment index](EXPERIMENTS.md) · [中文](EXPERIMENTS.zh-CN.md) | Stage-by-stage code and record navigation |
| [`simulator/cpp/`](simulator/cpp/) | Current C++20 main project: energies, implicit Euler, Projected Newton, meshes and scenes, executable examples, and small test inputs |
| [`simulator/Python/`](simulator/Python/) | Historical Python reference implementation; its former import name is described in the included notes |
| [`simulator/tools/`](simulator/tools/) | Structured-rod mesh generation, S2 runner, plotting, and related tools |
| [`experiments/workspaces/`](experiments/workspaces/) | Experiment source and scripts separated by run; each directory retains its own implementation without replacing the main project |
| [`simulator/docs/`](simulator/docs/) | Specifications, reports, and diagnostic records for S0/T0 through S6 |
| [`planning/`](planning/) | Historical validation protocols |
| [`assets/ipc/`](assets/ipc/) | Small torus and sphere1K meshes and their source records; other large meshes must be obtained separately |
| [`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json) | Original locations and SHA-256 hashes for files imported into this snapshot |

This snapshot was assembled from a local `fiction` workspace. The production project’s original Git HEAD was `b3745128119eb5def111159dbd3310a9e31fcc6f`. The snapshot also includes working files that were uncommitted or untracked at the time, so it is not an exact checkout of that commit. The `cpp/` trees in different workspaces are isolated versions from different stages; reproduce an original experiment using the source snapshot, parameters, and inputs from its corresponding run.

## Build the main project

The project requires C++20, CMake 3.20+, and Eigen 3.4+; SuiteSparse/CHOLMOD is optional. On a machine without a graphics environment, the viewer can be disabled:

```bash
cmake -S simulator/cpp -B build -DCMAKE_BUILD_TYPE=Release -DFICTION_USE_VIEWER=OFF
cmake --build build -j 3
ctest --test-dir build --output-on-failure
```

Viewer, PNG, and video features use GLFW/OpenGL, zlib, and ffmpeg according to the build options. See the [experiment index](EXPERIMENTS.md) for stage-specific runners, commands, and design notes. Scripts are imported as plain text; shell scripts can be run with `bash path/to/script.sh`. Some historical scripts retain absolute paths and build locations from the original machine. Check these against your environment before use, and do not overwrite the original records.

## Inputs and third-party content

Generate structured square rods with [`simulator/tools/make_structured_bar.py`](simulator/tools/make_structured_bar.py). The torus and sphere1K meshes come from the IPC project; their original sources and hashes are recorded in [`assets/ipc/sources.json`](assets/ipc/sources.json). The rod meshes in `simulator/cpp/tests/data/` were converted from ElastoGen examples. The Armadillo13K, 28K, and 54K meshes are not included in this source snapshot; their file hashes and historical local input locations are listed in the [experiment index](EXPERIMENTS.md).

Source code exploring contact in isolated experiments is retained in the corresponding workspaces. The external IPC Toolkit, Catch2 source, and build caches are not included. Full original trajectories, logs, executables, and frozen source archives remain in the local fiction project. Selected article figures and videos are now included under `blog/`. Historical paths and run IDs in `simulator/docs/` are preserved as experiment records; some links require access to the original workspace.

No project source-code license has been specified. Some meshes come from third parties; [`assets/ipc/UPSTREAM_LICENSE`](assets/ipc/UPSTREAM_LICENSE) is a copy of the upstream repository's license. The project license and mesh permissions have not been resolved here.

## Check the article files

```bash
python tools/check_blog.py
```

This checks bilingual emphasis, mathematical expressions, figure numbering, paragraph alignment, relative links, and bundled file hashes. It does not run a simulation.
