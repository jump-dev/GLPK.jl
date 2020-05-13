# TODO(odow):
#
# This script can be used to build the C interface to Clp. However, it requires
# you to manually do the following steps first:
#
# 1) Copy glpk.h from GLPK into this /scripts directory
#
# It should be possible to build the wrapper using the jll's, but I couldn't
# figure out how to do that, and I didn't have enough time to spend on it.

import Clang

LIBCLP_HEADERS = [
    joinpath(@__DIR__, "glpk.h"),
]

const COMMON = joinpath(dirname(@__DIR__), "src", "gen", "libglpk_common.jl")

wc = Clang.init(
    headers = LIBCLP_HEADERS,
    output_file = joinpath(dirname(@__DIR__), "src", "gen", "libglpk_api.jl"),
    common_file = COMMON,
    header_wrapped = (root, current) -> root == current,
    header_library = x -> "libglpk",
    clang_diagnostics = true,
)

run(wc)

function manual_corrections()
    file = read(COMMON, String)
    for correction in [
        "tail::Ptr{glp_vertex}" => "tail::Ptr",
        "head::Ptr{glp_vertex}" => "head::Ptr",
        "struct " => "mutable struct ",
    ]
        file = replace(file, correction)
    end
    for s_type in ["glp_smcp", "glp_iptcp", "glp_iocp"]
        r = Regex("(mutable struct $(s_type).+?end)", "s")
        str = replace(
            match(r, file)[1], "end" => "    $(s_type)() = new()\nend"
        )
        file = replace(file, r => str)
    end
    write(COMMON, file)
end
manual_corrections()

rm(joinpath(dirname(@__DIR__), "src", "gen", "LibTemplate.jl"))
