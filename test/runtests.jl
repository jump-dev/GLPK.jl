# Copyright (c) 2012 Carlo Baldassi and GLPK.jl contributors
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

using GLPK
using Test

@testset "C API" begin
    @testset "glpk_tst_$(i).jl" for i in 1:6
        include("glpk_tst_$i.jl")
    end
end

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
    include("MOI_callbacks.jl")
end
