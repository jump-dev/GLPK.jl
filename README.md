#  GLPK.jl

[![Build Status](https://github.com/jump-dev/GLPK.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/jump-dev/GLPK.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/GLPK.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/GLPK.jl)

[GLPK.jl](https://github.com/jump-dev/GLPK.jl) is a wrapper for the [GNU Linear Programming Kit library](https://www.gnu.org/software/glpk).

The wrapper has two components:

 * a thin wrapper around the complete C API
 * an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Affiliation

This wrapper is maintained by the JuMP community and is not an GNU project.

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please [open a GitHub issue](https://github.com/jump-dev/GLPK.jl/issues/new).

## License

`GLPK.jl` is licensed under the [GPL v3 license](https://github.com/jump-dev/GLPK.jl/blob/master/LICENSE.md).

## Installation

Install GLPK using `Pkg.add`:
```julia
import Pkg
Pkg.add("GLPK")
```

In addition to installing the GLPK.jl package, this will also download and
install the GLPK binaries. You do not need to install GLPK separately.

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Use with JuMP

To use GLPK with JuMP, use `GLPK.Optimizer`:
```julia
using JuMP, GLPK
model = Model(GLPK.Optimizer)
set_attribute(model, "tm_lim", 60 * 1_000)
set_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
```

If the model is primal or dual infeasible, GLPK will attempt to find a
certificate of infeasibility. This can be expensive, particularly if you do not
intend to use the certificate. If this is the case, use:
```julia
model = Model(() -> GLPK.Optimizer(; want_infeasibility_certificates = false))
```

## MathOptInterface API

The GLPK optimizer supports the following constraints and attributes.

List of supported objective functions:

 * [`MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}`](@ref)

List of supported variable types:

 * [`MOI.Reals`](@ref)

List of supported constraint types:

 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Integer`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.ZeroOne`](@ref)

List of supported model attributes:

 * [`MOI.HeuristicCallback()`](@ref)
 * [`MOI.LazyConstraintCallback()`](@ref)
 * [`MOI.Name()`](@ref)
 * [`MOI.ObjectiveSense()`](@ref)
 * [`MOI.UserCutCallback()`](@ref)

## Options

Options for GLPK are comprehensively documented in the [PDF documentation](https://github.com/jump-dev/GLPK.jl/files/11143880/glpk.pdf),
but they are hard to find.

 * Options when solving a linear program are defined in Section 2.8.1
 * Options when solving a mixed-integer program are defined in Section 2.10.5

However, the following options are likely to be the most useful:

| Parameter      | Example            | Explanation                            |
| -------------- | ------------------ | -------------------------------------- |
| `msg_lev`      | `GLPK.GLP_MSG_ALL` | Message level for terminal output      |
| `presolve`     | `GLPK.GLP_ON`      | Turn presolve on or off                |
| `tol_int`      | `1e-5`             | Absolute tolerance for integer feasibility |
| `tol_obj`      | `1e-7`             | Relative objective tolerance for mixed-integer programs |

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

## C API

The C API can be accessed via `GLPK.glp_XXX` functions, where the names and
arguments are identical to the C API. See the `/tests` folder for inspiration.

## Thread safety

GLPK is not thread-safe and should not be used with multithreading.
