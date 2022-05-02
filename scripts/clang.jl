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

import Clang
import GLPK_jll

GEN_DIR = joinpath(dirname(@__DIR__), "src", "gen")

const COMMON = joinpath(GEN_DIR, "libglpk_common.jl")

wc = Clang.init(
    headers = [
        joinpath(GLPK_jll.artifact_dir, "include", "glpk.h"),
    ],
    output_file = joinpath(GEN_DIR, "libglpk_api.jl"),
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

    function new_functions(name, constr)
        s = constr ? "\n    $(name)() = new()\nend" : "\nend"
        return """
        $(s)

        Base.cconvert(::Type{Ptr{$(name)}}, x::$(name)) = x
        function Base.unsafe_convert(::Type{Ptr{$(name)}}, x::$(name))
            return convert(Ptr{$(name)}, pointer_from_objref(x))
        end"""
    end
    for (s_type, constr) in [
        "glp_bfcp" => false,
        "glp_smcp" => true,
        "glp_iptcp" => true,
        "glp_iocp" => true,
        "glp_attr" => false,
        "glp_mpscp" => true,
        "glp_cpxcp" => true,
        "glp_arc" => false,
        "glp_vertex" => false,
        "glp_graph" => false,
    ]
        r = Regex("(mutable struct $(s_type).+?end)", "s")
        str = replace(match(r, file)[1], "\nend" => new_functions(s_type, constr))
        file = replace(file, r => str)
    end
    write(COMMON, file)
end
manual_corrections()

rm(joinpath(GEN_DIR, "LibTemplate.jl"))
