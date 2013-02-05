module GLPK_tests

function glpk_tst_all()
    prev_len = 0
    for i = 1:6
        f = "glpk_tst_$i.jl"
        println("Running $f")
        include(f)
    end
end

end

GLPK_tests.glpk_tst_all()
