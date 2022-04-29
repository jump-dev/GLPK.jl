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

using Test
import GLPK

# Testing reading and writing problem files

function glpk_tst_2()
    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    lp = GLPK.glp_create_prob()

    # test GLPK native format
    GLPK.glp_read_prob(lp, 0, joinpath(datadir, "sample.prob"))
    @test GLPK.glp_simplex(lp, C_NULL) == 0

    @test GLPK.glp_read_prob(lp, 0, joinpath(datadir, "plan.lp")) != 0
    @test GLPK.glp_write_prob(lp, 0, "") != 0

    mktemp() do path, _
        return GLPK.glp_write_prob(lp, 0, path)
    end

    # test MPS format
    @test GLPK.glp_read_mps(
        lp,
        GLPK.GLP_MPS_FILE,
        C_NULL,
        "nonexisting_file",
    ) != 0

    GLPK.glp_read_mps(
        lp,
        GLPK.GLP_MPS_DECK,
        C_NULL,
        joinpath(datadir, "plan.mps"),
    )
    @test GLPK.glp_simplex(lp, C_NULL) == 0

    @test GLPK.glp_read_mps(
        lp,
        GLPK.GLP_MPS_DECK,
        C_NULL,
        joinpath(datadir, "plan.lp"),
    ) != 0
    @test GLPK.glp_write_mps(lp, GLPK.GLP_MPS_FILE, C_NULL, "") != 0

    mktemp() do path, _
        return GLPK.glp_write_mps(lp, GLPK.GLP_MPS_FILE, C_NULL, path)
    end

    # Test LP format
    @test GLPK.glp_read_lp(lp, C_NULL, "nonexisting_file") != 0

    GLPK.glp_read_lp(lp, C_NULL, joinpath(datadir, "plan.lp"))
    @test GLPK.glp_simplex(lp, C_NULL) == 0

    @test GLPK.glp_read_lp(lp, C_NULL, joinpath(datadir, "plan.mps")) != 0
    @test GLPK.glp_write_lp(lp, C_NULL, "") != 0

    mktemp() do path, _
        return GLPK.glp_write_lp(lp, C_NULL, path)
    end

    GLPK.glp_term_out(prev_term_out)
    return GLPK.glp_delete_prob(lp)
end

glpk_tst_2()
