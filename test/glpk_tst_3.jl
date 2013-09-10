using Base.Test
import GLPK

# Test reading and writing MathProg files

function glpk_tst_3()
    prev_term_out = GLPK.term_out(GLPK.OFF)

    datadir = joinpath(Pkg.dir(), "GLPK", "test", "data")
    isdir(datadir) || (datadir = joinpath(Pkg.dir(), "GLPK.jl", "test", "data"))

    lp = GLPK.Prob()
    tran = GLPK.MathProgWorkspace()

    GLPK.mpl_read_model(tran, joinpath(datadir, "queens.mod"), 0)

    outfile = tempname()
    try
        GLPK.mpl_generate(tran, outfile)
    finally
        isfile(outfile) && rm(outfile)
    end
    GLPK.mpl_build_prob(tran, lp)

    @test GLPK.simplex(lp) == 0

    GLPK.term_out(prev_term_out)
end

glpk_tst_3()
