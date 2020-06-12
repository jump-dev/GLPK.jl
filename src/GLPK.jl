module GLPK

if haskey(ENV,"JULIA_GLPK_LIBRARY_PATH") || VERSION < v"1.3"
    deps_file = joinpath(dirname(@__DIR__), "deps", "deps.jl")
    if isfile(deps_file)
        include(deps_file)
    else
        error("GLPK not properly installed. Please run import `Pkg; Pkg.build(\"GLPK\")`.")
    end
else
    import GLPK_jll: libglpk
end

using CEnum

include("gen/ctypes.jl")
include("gen/libglpk_common.jl")
include("gen/libglpk_api.jl")

const _GLPK_VERSION = VersionNumber("$GLP_MAJOR_VERSION.$GLP_MINOR_VERSION.0")

if !(v"4.64.0" <= _GLPK_VERSION <= v"4.64.0")
    error(
        "You have installed version $_GLPK_VERSION of GLPK, which is not " *
        "supported by GLPK.jl. If the version change was breaking, changes " *
        "will need to be made to the Julia code. Please open an issue at " *
        "https://github.com/JuliaOpt/GLPK.jl."
    )
end

include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_callbacks.jl")

# GLPK exports everything except internal symbols, which are defined as those
# whose name starts with an underscore. If you don't want all of these symbols
# in your environment, then use `import GLPK` instead of `using GLPK`.

# Do not add GLPK-defined symbols to this exclude list. Instead, rename them
# with an underscore.
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all=true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_")
        continue
    end
    if !(Base.isidentifier(sym) || (startswith(sym_string, "@") &&
         Base.isidentifier(sym_string[2:end])))
       continue
    end
    @eval export $sym
end

end
