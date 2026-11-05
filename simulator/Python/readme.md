# Mass-Spring Solids Simulation

> Historical reference only: the files still import the old package name
> `Python_slovement`, while this directory is named `Python`. The commands below
> describe the original layout and do not run unchanged in this checkout. For
> the supported build, tests, and self-contained example, use the
> [repository README](../../README.md#快速开始). No Python dependency installation
> is needed for the C++ quick start.

A square falling onto a ground under gravity and then compressed by a ceiling is simulated with mass-spring elasticity potential and implicit Euler time integration.
Each time step is solved by minimizing the Incremental Potential with the projected Newton method.

## Dependencies
```
pip install numpy scipy pygame
```

## Run
```
python simulator.py
```
