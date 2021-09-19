module GLPK

if haskey(ENV, "JULIA_GLPK_LIBRARY_PATH") || VERSION < v"1.3"
    deps_file = joinpath(dirname(@__DIR__), "deps", "deps.jl")
    if isfile(deps_file)
        include(deps_file)
    else
        error(
            "GLPK not properly installed. Please run " *
            "`import Pkg; Pkg.build(\"GLPK\")`.",
        )
    end
else
    import GLPK_jll: libglpk
end

using CEnum

include("gen/ctypes.jl")
include("gen/libglpk_common.jl")
include("gen/libglpk_api.jl")

const GLP_DBL_MAX = prevfloat(Inf)

"""
    offset(x::Vector)

GLPK uses 1-based indexing for its arrays. But since C has 0-based indexing, all
1-based vectors passed to GLPK need to be padded with a "0'th" element that will
never be accessed. To avoid doing this padding in Julia, we convert the vector
to a reference, and use the optional second argument to ensure the reference
points to the "0'th" element of the array. This is safe to do, provided C never
accesses `x[0]`.

See the GLPK manual for more details.
"""
offset(x::Vector) = Ref(x, 0)

const _GLPK_VERSION = let
    p = glp_version()
    VersionNumber(parse.(Int, split(unsafe_string(p), "."))...)
end

if !(v"4.64.0" <= _GLPK_VERSION <= v"5.0.0")
    error("""
    You have installed version $_GLPK_VERSION of GLPK, which is not supported
    by GLPK.jl. We prefer GLPK version 5.0, but this may also work on versions
    4.64 and newer.

    After installing GLPK 5.0, run:

        import Pkg
        Pkg.rm("GLPK")
        Pkg.add("GLPK")

    If you have a newer version of GLPK installed, changes may need to be made
    to the Julia code. Please open an issue at
    https://github.com/jump-dev/GLPK.jl.
    """)
end

include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_copy.jl")
include("MOI_wrapper/MOI_callbacks.jl")
include("MOI_wrapper/deprecated_constants.jl")

# GLPK exports all `GLP_XXX` and `glp_xxx` symbols. If you don't want all of
# these symbols in your environment, then use `import GLPK` instead of
# `using GLPK`.

for sym in names(@__MODULE__, all = true)
    sym_string = string(sym)
    if startswith(sym_string, "GLP_") || startswith(sym_string, "glp_")
        @eval export $sym
    end
end

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end
