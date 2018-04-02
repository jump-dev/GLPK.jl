using Compat.Test
import GLPK

# Testing reading and writing problem files

function glpk_tst_2()
    prev_term_out = GLPK.term_out(GLPK.OFF)

    datadir = joinpath(dirname(@__FILE__), "data")
    @assert isdir(datadir)

    lp = GLPK.Prob()

    # test GLPK native format
    GLPK.read_prob(lp, 0, joinpath(datadir, "sample.prob"))
    @test GLPK.simplex(lp) == 0

    @glpk_test_throws GLPK.read_prob(lp, 0, joinpath(datadir, "plan.lp"))
    @glpk_test_throws GLPK.write_prob(lp, 0, "") != 0

    filecopy = tempname()

    try
        GLPK.write_prob(lp, 0, filecopy)
    finally
        isfile(filecopy) && rm(filecopy)
    end

    # test MPS format
    @glpk_test_throws GLPK.read_mps(lp, GLPK.MPS_FILE, nothing, "nonexisting_file")
    @glpk_test_throws GLPK.read_mps(lp, GLPK.MPS_FILE, C_NULL, "nonexisting_file")
    @glpk_test_throws GLPK.read_mps(lp, GLPK.MPS_FILE, "nonexisting_file")

    GLPK.read_mps(lp, GLPK.MPS_DECK, C_NULL, joinpath(datadir, "plan.mps"))
    @test GLPK.simplex(lp, nothing) == 0

    @glpk_test_throws GLPK.read_mps(lp, GLPK.MPS_DECK, C_NULL, joinpath(datadir, "plan.lp"))
    @glpk_test_throws GLPK.write_mps(lp, GLPK.MPS_FILE, C_NULL, "")
    @glpk_test_throws GLPK.write_mps(lp, GLPK.MPS_FILE, nothing, "")
    @glpk_test_throws GLPK.write_mps(lp, GLPK.MPS_FILE, "")

    try
        GLPK.write_mps(lp, GLPK.MPS_FILE, C_NULL, filecopy)
    finally
        isfile(filecopy) && rm(filecopy)
    end

    # Test LP format
    @glpk_test_throws GLPK.read_lp(lp, C_NULL, "nonexisting_file")

    GLPK.read_lp(lp, C_NULL, joinpath(datadir, "plan.lp"))
    @test GLPK.simplex(lp) == 0

    @glpk_test_throws GLPK.read_lp(lp, C_NULL, joinpath(datadir, "plan.mps"))
    @glpk_test_throws GLPK.write_lp(lp, C_NULL, "")

    try
        GLPK.write_lp(lp, C_NULL, filecopy)
    finally
        isfile(filecopy) && rm(filecopy)
    end
    
    GLPK.term_out(prev_term_out)
end

glpk_tst_2()
