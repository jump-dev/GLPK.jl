function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(MathOptInterface.optimize!),Optimizer})   # time: 0.17121503
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.ScalarAffineFunction{Float64},Union{MathOptInterface.EqualTo{Float64}, MathOptInterface.GreaterThan{Float64}, MathOptInterface.LessThan{Float64}}})   # time: 0.072197944
    Base.precompile(Tuple{typeof(MathOptInterface.add_variables),Optimizer,Int64})   # time: 0.049029674
    Base.precompile(Tuple{typeof(MathOptInterface.submit),Optimizer,MathOptInterface.HeuristicSolution{CallbackData},Vector{MathOptInterface.VariableIndex},Vector{Float64}})   # time: 0.017970324
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,MathOptInterface.Silent,Bool})   # time: 0.015073183
    Base.precompile(Tuple{typeof(MathOptInterface.modify),Optimizer,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}},MathOptInterface.ScalarCoefficientChange{Float64}})   # time: 0.010970611
    Base.precompile(Tuple{typeof(MathOptInterface.delete),Optimizer,MathOptInterface.VariableIndex})   # time: 0.010189695
    Base.precompile(Tuple{typeof(MathOptInterface.submit),Optimizer,MathOptInterface.LazyConstraint{CallbackData},MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}})   # time: 0.009361999
    # TODO: Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,MathOptInterface.ObjectiveFunction{F<:MathOptInterface.ScalarAffineFunction{Core.Float64}},F<:MathOptInterface.ScalarAffineFunction{Core.Float64}})   # time: 0.00935949
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.VariablePrimal,MathOptInterface.VariableIndex})   # time: 0.008826303
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}})   # time: 0.007045114
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.SingleVariable,MathOptInterface.EqualTo{Float64}})   # time: 0.005506918
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.SingleVariable,MathOptInterface.LessThan{Float64}})   # time: 0.005424305
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.ConstraintPrimal,Union{MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{T}} where T, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{T}} where T}})   # time: 0.005108291
    Base.precompile(Tuple{typeof(MathOptInterface.submit),Optimizer,MathOptInterface.UserCut{CallbackData},MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}})   # time: 0.004174831
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.ConstraintDual,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}})   # time: 0.00381129
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,CallbackFunction,Function})   # time: 0.003472098
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.SingleVariable,MathOptInterface.Integer})   # time: 0.003090693
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.TerminationStatus})   # time: 0.002651099
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.CallbackNodeStatus{CallbackData}})   # time: 0.002032824
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,MathOptInterface.ObjectiveFunction{MathOptInterface.ScalarAffineFunction{Float64}},MathOptInterface.ScalarAffineFunction{Float64}})   # time: 0.002016101
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.ObjectiveValue})   # time: 0.001891955
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.PrimalStatus})   # time: 0.001678843
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}})   # time: 0.001646429
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}})   # time: 0.001606696
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.DualStatus})   # time: 0.001574324
    Base.precompile(Tuple{typeof(MathOptInterface.is_valid),Optimizer,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable, MathOptInterface.ZeroOne}})   # time: 0.001551993
    Base.precompile(Tuple{typeof(MathOptInterface.add_constraint),Optimizer,MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}})   # time: 0.001549463
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,MathOptInterface.LazyConstraintCallback,Function})   # time: 0.001537747
    Base.precompile(Tuple{typeof(MathOptInterface.is_valid),Optimizer,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable, MathOptInterface.Integer}})   # time: 0.001446975
    Base.precompile(Tuple{typeof(MathOptInterface.add_variable),Optimizer})   # time: 0.001420942
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,MathOptInterface.ObjectiveSense,MathOptInterface.OptimizationSense})   # time: 0.001282038
end
