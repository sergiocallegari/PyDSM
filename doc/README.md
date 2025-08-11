# Building the pydsm documetation

(Re)generating the documentation requires:

- `sphinx`
- `numpydoc`
    - extension used by `sphinx` to manage docstrings according to the NumPy conventions

Build the documentation from this directory using the makefile. E.g., `make html` or more probably something like `uv run make html`.
