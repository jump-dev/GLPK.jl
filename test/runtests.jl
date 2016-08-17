module GLPK_tests

macro test_throws_02(args...)
    :(@test_throws($(esc(args[1])), $(esc(args[2]))))
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
