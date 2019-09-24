using GLPK, Test

macro test_throws_02(args...)
    :(@test_throws($(esc(args[1])), $(esc(args[2]))))
end

macro glpk_test_throws(args)
    :(@test_throws_02 GLPK.GLPKError $(esc(args)))
end

@testset "C API" begin
    for i = 1:6
        f = "glpk_tst_$i.jl"
        @testset "$f" begin
            include(f)
        end
    end
end

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
    include("MOI_callbacks.jl")
end
