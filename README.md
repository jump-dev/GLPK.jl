#  GLPK.jl

| **Build Status**                                                                                    |
|:---------------------------------------------------------------------------------------------------:|
| [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coveralls-img]][coveralls-url] |


GLPK.jl is a wrapper for the [GNU Linear Programming Kit library](http://www.gnu.org/software/glpk).
It makes it possible to access nearly all of GLPK functionality from within Julia programs.

This package is part of [the JuliaOpt project](http://www.juliaopt.org/).

## Installation

The package is registered in the [General registry](https://github.com/JuliaRegistries/General/) and so can be installed with `Pkg.add`.

```julia
julia> import Pkg

julia> Pkg.add("GLPK")
```

GLPK.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl) to automatically install the GLPK binaries with [GMP](https://gmplib.org) support.

## Custom Installation

After GLPK.jl is installed and built, you can replace the installed binary dependencies with custom builds by overwriting the binaries and libraries in GLPK.jl's `deps/usr` folder. For instance, Julia v1.0 and GLPK v0.9.1 this can be achieved by running
```bash
./configure --prefix=$(julia -e 'import GLPK; println(dirname(dirname(pathof(GLPK))))')/deps/usr
make
make install
```
in GLPK's source folder.

Note that the custom binaries will not be overwritten by subsequent builds of the currently installed version of GLPK.jl. However, if GLPK.jl is updated and the update includes new BinaryProvider versions of the GLPK binaries, then the custom binaries will be overwritten by the new BinaryProvider versions.

## `GLPK.Optimizer`

Use `GLPK.Optimizer` to create a new optimizer object:
```julia
using GLPK
model = GLPK.Optimizer(tm_lim = 60.0, msg_lev = GLPK.OFF)
```
For JuMP, use:
```julia
using JuMP, GLPK
model = Model(
    with_optimizer(GLPK.Optimizer, tm_lim = 60.0, msg_lev = GLPK.OFF)
)
```

**Note: previous versions of `GLPK.jl` required you to choose either `GLPKSolverLP` or `GLPKSolverMIP`. This is no longer needed; just use `GLPK.Optimizer`.**

## Pre-emptive checks

`GLPK.jl` has a lot of pre-emptive checks to catch cases where the C API might
throw an uninformative error. However, in benchmarks, this takes a
non-negligible amount of time (e.g. 20% in add_constraints). At the risk of
possibly running into an uninformative error, you can run the following after
importing GLPK to disable these checks:
```julia
using GLPK
GLPK.jl_set_preemptive_check(false)
```

[travis-img]: https://api.travis-ci.org/JuliaOpt/GLPK.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaOpt/GLPK.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/4t5e2dir3gp7fb6h?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaOpt/glpk-jl

[coveralls-img]: https://img.shields.io/coveralls/JuliaOpt/GLPK.jl.svg
[coveralls-url]: https://coveralls.io/r/JuliaOpt/GLPK.jl
