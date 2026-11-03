# Article and figure sources

The [Chinese](fiction.md) and [English](fiction-en.md) articles share the same 21 numbered figures, three saved videos, and one animation that switches between prescribed loading cases. The articles remain drafts for author review.

## How the figures were made

| Figures | Content and source | Type |
|---|---|---|
| 1 | A particle released from rest; gravity 9.81 m/s² and a 0.005 s time step give continuous and implicit-Euler displacements of 0.122625 and 0.245250 mm | Analytic illustration |
| 2 | Reference, stretch, shear, and compression configurations in [`mesh_check.cpp`](../simulator/cpp/examples/mesh_check.cpp) | Prescribed element geometry |
| 3–4 | Linear oscillator scene, phase trajectories, state differences, and energy retention; `stage0-20261030-01` | Saved trajectory and scene schematic |
| 5–6 | Uniform rod stretch; lateral reference factors and reactions from `s1-merged-20261031-01` | Prescribed affine reference geometry and saved result summaries |
| 7–8 | Static cantilever shapes and mesh-refinement data; `s2-extension-20261031-01` | Saved node states and result summaries |
| 9–12 | Dynamic cantilever motion, saved video, and time-step trajectories; `s3-visual-main-20261102-01`, `s3-bend-20261031-01`, `s3-dtmesh-20261102-01`, `s3-dtmesh-ref-20261102-01` | Saved geometry, frames, and trajectories |
| 13–15 | Torsion geometry, cross-sections, video, rotation, and torque; `s3t-twist-20261031-01` and `s3t-visual-main-20261102-02` | Saved states and trajectories |
| 16–17 | Torus geometry at 0, 1, and 2 s and the saved video; `s4-nonrod-20261031-01` | Saved surface meshes and frames |
| 18–19 | Armadillo reference geometry, fixed-node mask, timing, and volume ratios; `s5-scale-20261031-01` and `s5-54k-main-20261102-01` | Reference geometry and saved run measurements |
| 20–21 | Four-body translation layout, block structure, and sixteen-body states; `s6a-four-main-20261102-01`, `s6a-sixteen-main-20261102-01`, `s6a-layout-main-20261102-01` | Saved positions and a matrix-structure schematic |

Figure 13 enlarges transverse dimensions by 3× in its left panel; its cross-section panel uses actual coordinates. The loading-case GIF for Figure 5 shows prescribed reference shapes, not physical time integration. Figure 18 uses native Armadillo coordinates, with +y upright and gravity along −z. Figure 20 uses circles as position symbols and illustrates matrix blocks rather than actual sparse entries.

PNG files and videos are in [`assets/`](assets/). The ten new scene illustrations also have SVG originals in [`editable-figures/`](editable-figures/). The complex Armadillo surface is rasterized inside its SVG; labels and annotations remain vector elements.

## Reports and traceability

- [Material and time integration](../simulator/docs/stage0-validation-20261030.md)
- [Uniform stretch](../simulator/docs/s1-reacceptance-20261031.md)
- [Static cantilever refinement](../simulator/docs/s2-extension-results-20261031.md)
- [Dynamic cantilever](../simulator/docs/s3-bend-results-20261031.md) and [four-time-step addendum](../simulator/docs/s3-dtmesh-000625-addendum-20261102.md)
- [Torsion](../simulator/docs/s3t-twist-results-20261031.md), [torus](../simulator/docs/s4-nonrod-results-20261031.md), and [54K scale run](../simulator/docs/s5-54k-main-results-20261102.md)
- [Four bodies](../simulator/docs/s6a-four-main-results-20261102.md) and [sixteen bodies](../simulator/docs/s6a-sixteen-main-results-20261102.md)
- [Stopping-criterion integration](../simulator/docs/stopcriterion-merge-20261102.md)

The original CSV files, full node states, frame sequences, and frozen run archives remain in the local research archive; they are not all bundled in this repository. Run IDs above identify those inputs. [`figure-provenance.json`](figure-provenance.json) retains input identifiers and hashes, using archive references rather than machine-specific absolute paths. [`MANIFEST.sha256`](MANIFEST.sha256) covers the files distributed with the article.

## Formatting and language

Chinese prose paragraphs begin with two ideographic spaces. English paragraphs are flush left. The 18 emphasized spans correspond between languages, and all 168 mathematical expressions are preserved. Both versions use the same translation key and remain `draft: true`.

The source code and experiment archive retain their internal stage IDs for traceability. The article uses descriptive experiment titles throughout.
