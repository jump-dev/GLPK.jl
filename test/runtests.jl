module GLPK_tests

# backwards-compatible test_throws (works in julia 0.2)
macro test_throws_02(args...)
    if VERSION >= v"0.3-"
        :(@test_throws($(esc(args[1])), $(esc(args[2]))))
    else
        :(@test_throws($(esc(args[2]))))
    end
end

macro glpk_test_throws(args)
    :(@test_throws_02 GLPK.GLPKError $(esc(args)))
end

for i = 1:6
    f = "glpk_tst_$i.jl"
    println("Running $f")
    include(f)
end

end
