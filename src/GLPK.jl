# Copyright (c) 2012 GLPK.jl contributors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the Licence, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

module GLPK

import GLPK_jll

function __init__()
    global libglpk = GLPK_jll.libglpk
    p = glp_version()
    version = VersionNumber(parse.(Int, split(unsafe_string(p), "."))...)
    if !(v"4.64.0" <= version <= v"5.0.0")
        error("""
        You have installed version $version of GLPK, which is not supported
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
    return
end

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

include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_copy.jl")
include("MOI_wrapper/MOI_callbacks.jl")

# GLPK exports all `GLP_XXX` and `glp_xxx` symbols. If you don't want all of
# these symbols in your environment, then use `import GLPK` instead of
# `using GLPK`.

for sym in names(@__MODULE__, all = true)
    sym_string = string(sym)
    if startswith(sym_string, "GLP_") || startswith(sym_string, "glp_")
        @eval export $sym
    end
end

include("precompile.jl")
_precompile_()

end
