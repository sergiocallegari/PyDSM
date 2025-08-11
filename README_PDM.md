# Use of `pdm` as a project management tool

The package is developed using [uv](https://github.com/astral-sh/uv) for project management.

Building the project with [PDM](https://github.com/pdm-project/pdm) should also work, but is not very much tested. Furthermore, some caution is required for the build. Specifically, the required versions of `scs` and `cvxpy` must currently be added using the `--no-isolation` flag.
