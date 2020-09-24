#  GLPK.jl

| **Build Status**                                                                                    |
|:---------------------------------------------------------------------------------------------------:|
| [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coveralls-img]][coveralls-url] |

GLPK.jl is a wrapper for the [GNU Linear Programming Kit library](https://www.gnu.org/software/glpk).

It has two components:
 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

The C API can be accessed via `GLPK.glp_XXX` functions, where the names and
arguments are identical to the C API. See the `/tests` folder for inspiration.

## Installation

The package is registered in the [General registry](https://github.com/JuliaRegistries/General/)
and so can be installed with `Pkg.add`.

```julia
import Pkg
Pkg.add("GLPK")
```

In addition to installing the GLPK.jl package, this will also download and
install the GLPK binaries. (You do not need to install GLPK separately.) If you
require a custom build of GLPK, see the **Custom Installation** instructions
below.

## Custom Installation

To install custom built GLPK binaries, use the environmental variable
`JULIA_GLPK_LIBRARY_PATH` to point to the path at which you installed GLPK (the
folder containing `libglpk`). For example, on Mac, after installing GLPK with
`brew install glpk`, use:
```julia
ENV["JULIA_GLPK_LIBRARY_PATH"] = "/usr/local/Cellar/glpk/4.65/lib"
import Pkg
Pkg.add("GLPK")
Pkg.build("GLPK")
```
Replace `"/usr/local/Cellar/glpk/4.65/lib"` with a different path as
appropriate.

**You must have `JULIA_GLPK_LIBRARY_PATH` set _every_ time you run `using GLPK`,
not just when you install it.**

Switch back to the default binaries as follows:
```julia
delete!(ENV, "JULIA_GLPK_LIBRARY_PATH")
import Pkg
Pkg.build("GLPK")
```

## Use with JuMP

We highly recommend that you use GLPK.jl with higher level packages such as
[JuMP.jl](https://github.com/jump-dev/JuMP.jl).

This can be done using the `GLPK.Optimizer` object. Here is how to create a JuMP
model that uses GLPK as the solver.

```julia
using JuMP, GLPK

model = Model(GLPK.Optimizer)
set_optimizer_attribute(model, "tm_lim", 60 * 1_000)
set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
```

## Callbacks

Here is an example using GLPK's solver-specific callbacks.

```julia
using JuMP, GLPK, Test

model = Model(GLPK.Optimizer)
@variable(model, 0 <= x <= 2.5, Int)
@variable(model, 0 <= y <= 2.5, Int)
@objective(model, Max, y)
reasons = UInt8[]
function my_callback_function(cb_data)
    reason = GLPK.glp_ios_reason(cb_data.tree)
    push!(reasons, reason)
    if reason != GLPK.GLP_IROWGEN
        return
    end
    x_val = callback_value(cb_data, x)
    y_val = callback_value(cb_data, y)
    if y_val - x_val > 1 + 1e-6
        con = @build_constraint(y - x <= 1)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    elseif y_val + x_val > 3 + 1e-6
        con = @build_constraint(y - x <= 1)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
end
MOI.set(model, GLPK.CallbackFunction(), my_callback_function)
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL
@test primal_status(model) == MOI.FEASIBLE_POINT
@test value(x) == 1
@test value(y) == 2
@show reasons
```

[travis-img]: https://api.travis-ci.org/jump-dev/GLPK.jl.svg?branch=master
[travis-url]: https://travis-ci.org/jump-dev/GLPK.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/4t5e2dir3gp7fb6h?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaOpt/glpk-jl

[coveralls-img]: https://img.shields.io/coveralls/jump-dev/GLPK.jl.svg
[coveralls-url]: https://coveralls.io/r/jump-dev/GLPK.jl
