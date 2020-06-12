#  GLPK.jl

| **Build Status**                                                                                    |
|:---------------------------------------------------------------------------------------------------:|
| [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coveralls-img]][coveralls-url] |

GLPK.jl is a wrapper for the [GNU Linear Programming Kit library](http://www.gnu.org/software/glpk).

It has two components:
 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

The C API can be accessed via `GLPK.glp_XXX` functions, where the names and
arguments are identical to the C API. See the `/tests` folder for inspiration.

This package is part of [the JuliaOpt project](http://www.juliaopt.org/).

## Installation

The package is registered in the [General registry](https://github.com/JuliaRegistries/General/)
and so can be installed with `Pkg.add`.

```julia
julia> import Pkg

julia> Pkg.add("GLPK")
```

In addition to installing the GLPK.jl package, this will also download and
install the GLPK binaries. (You do not need to install GLPK separately.) If you
require a custom build of GLPK, see the instructions below.

## Custom Installation

To install custom built GLPK binaries set the environmental variable
`JULIA_GLPK_LIBRARY_PATH` and call `import Pkg; Pkg.build("GLPK")`. For example,
if the libraries are installed in `/opt/lib` just call:
```julia
ENV["JULIA_GLPK_LIBRARY_PATH"]="/opt/lib"
Pkg.build("GLPK")
```

If you do not want BinaryProvider to download the default binaries on install
set  `JULIA_GLPK_LIBRARY_PATH`  before calling `import Pkg; Pkg.add("GLPK")`.

To switch back to the default binaries clear `JULIA_GLPK_LIBRARY_PATH` and call
`import Pkg; Pkg.build("GLPK")`.

## `GLPK.Optimizer`

Use `GLPK.Optimizer` to create a new optimizer object for use with
MathOptInterface or JuMP.

In the following examples the time limit is set to one minute and logging is
turned off.
```julia
using GLPK
model = GLPK.Optimizer()
MOI.set(model, MOI.RawParameter("tm_lim"), 60_000) #  tm_lim is in milliseconds.
MOI.set(model, MOI.RawParameter("msg_lev"), GLPK.GLP_MSG_OFF)
```

For JuMP, use:
```julia
using JuMP, GLPK
model = Model(GLPK.Optimizer)
set_optimizer_attribute(model, "tm_lim", 60 * 1_000)
set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
```

**Note: previous versions of `GLPK.jl` required you to choose either `GLPKSolverLP` or `GLPKSolverMIP`. This is no longer needed; just use `GLPK.Optimizer`.**

[travis-img]: https://api.travis-ci.org/JuliaOpt/GLPK.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaOpt/GLPK.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/4t5e2dir3gp7fb6h?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaOpt/glpk-jl

[coveralls-img]: https://img.shields.io/coveralls/JuliaOpt/GLPK.jl.svg
[coveralls-url]: https://coveralls.io/r/JuliaOpt/GLPK.jl
