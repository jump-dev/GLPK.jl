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

using Test
import GLPK

# Test reading and writing MathProg files

function glpk_tst_3()
    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    lp = GLPK.glp_create_prob()
    tran = GLPK.glp_mpl_alloc_wksp()

    GLPK.glp_mpl_read_model(tran, joinpath(datadir, "queens.mod"), 0)

    outfile = tempname()
    try
        GLPK.glp_mpl_generate(tran, outfile)
        GLPK.glp_mpl_build_prob(tran, lp)
        GLPK.glp_mpl_free_wksp(tran)
    finally
        isfile(outfile) && rm(outfile)
    end

    @test GLPK.glp_simplex(lp, C_NULL) == 0

    GLPK.glp_term_out(prev_term_out)
    return GLPK.glp_delete_prob(lp)
end

glpk_tst_3()
