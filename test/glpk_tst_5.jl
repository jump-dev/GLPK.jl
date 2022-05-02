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

function glpk_tst_5()
    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    prev_term_out = GLPK.glp_term_out(GLPK.GLP_OFF)

    lp = GLPK.glp_create_prob()

    cnffile = joinpath(datadir, "sat_tst.cnf")
    GLPK.glp_read_cnfsat(lp, cnffile)
    GLPK.glp_check_cnfsat(lp)
    cnfcopy = joinpath(datadir, "sat_tst-copy.cnf")
    try
        GLPK.glp_write_cnfsat(lp, cnfcopy)
    finally
        isfile(cnfcopy) && rm(cnfcopy)
    end
    if Sys.WORD_SIZE != sizeof(Csize_t) * 8
        # changed condition to:
        # if (sizeof(void *) != sizeof(size_t))
        expected_ret = GLPK.GLP_EFAIL
    else
        expected_ret = 0
    end
    @test GLPK.glp_minisat1(lp) == expected_ret

    @test GLPK.glp_intfeas1(lp, 0, 0) == expected_ret
    return GLPK.glp_delete_prob(lp)
end

glpk_tst_5()
