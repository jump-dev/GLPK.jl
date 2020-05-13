using GLPK
using Test

@testset "C API" begin
    @testset "glpk_tst_$(i).jl" for i = 1:6
        if i ∈ [1, 2, 4]
            # TODO(odow): update these files.
            continue
        end
        include("glpk_tst_$i.jl")
    end
end

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
    include("MOI_callbacks.jl")
end
