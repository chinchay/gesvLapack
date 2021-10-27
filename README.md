# gesvLapack

`A * X = B` is solved by the `gesv` routine in LAPACK. Here, a Julia implementation is provided using several assumptions.

## Assumptions
- `A` and `B` are **SQUARE** matrices, bothe of the order `n`.

## Goals
- Reduce the number of allocations, and time solving

## Where to use it
- Intended for iterative processes, e.g. self-consistent equations. See tests directory.
