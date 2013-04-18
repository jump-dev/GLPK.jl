using Base.Test
import GLPK

# Testing reading and writing problem files

function glpk_tst_2()
    prev_term_out = GLPK.term_out(GLPK.OFF)

    datadir = joinpath(Pkg.dir(), "GLPK", "test", "data")

    lp = GLPK.Prob()

    # test GLPK native format
    GLPK.read_prob(lp, 0, joinpath(datadir, "sample.prob"))
    @test GLPK.simplex(lp) == 0

    @test_fails GLPK.read_prob(lp, 0, joinpath(datadir, "plan.lp"))
    @test_fails GLPK.write_prob(lp, 0, "") != 0

    filecopy = tempname()

    try
        GLPK.write_prob(lp, 0, filecopy)
    finally
        if isfile(filecopy)
            rm(filecopy)
        end
    end

    # test MPS format
    @test_fails GLPK.read_mps(lp, GLPK.MPS_FILE, C_NULL, "nonexisting_file")

    GLPK.read_mps(lp, GLPK.MPS_DECK, C_NULL, joinpath(datadir, "plan.mps"))
    @test GLPK.simplex(lp, nothing) == 0

    @test_fails GLPK.read_mps(lp, GLPK.MPS_DECK, C_NULL, joinpath(datadir, "plan.lp"))
    @test_fails GLPK.write_mps(lp, GLPK.MPS_FILE, C_NULL, "")

    try
        GLPK.write_mps(lp, GLPK.MPS_FILE, C_NULL, filecopy)
    finally
        if isfile(filecopy)
            rm(filecopy)
        end
    end

    # Test LP format
    @test_fails GLPK.read_lp(lp, C_NULL, "nonexisting_file")

    GLPK.read_lp(lp, C_NULL, joinpath(datadir, "plan.lp"))
    @test GLPK.simplex(lp) == 0

    @test_fails GLPK.read_lp(lp, C_NULL, joinpath(datadir, "plan.mps"))
    @test_fails GLPK.write_lp(lp, C_NULL, "")

    try
        GLPK.write_lp(lp, C_NULL, filecopy)
    finally
        if isfile(filecopy)
            rm(filecopy)
        end
    end
    
    GLPK.term_out(prev_term_out)
end

glpk_tst_2()
