function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    Base.precompile(MOI.set, (Optimizer, MOI.Silent, Bool))
    Base.precompile(MOI.copy_to, (Optimizer, MOI.Utilities.Model{Float64}))
    Base.precompile(MOI.add_variables, (Optimizer, Int))
    Base.precompile(MOI.add_variable, (Optimizer,))
    Base.precompile(MOI.delete, (Optimizer, MOI.VariableIndex))
    Base.precompile(MOI.get, (Optimizer, MOI.VariablePrimal, MOI.VariableIndex))

    Base.precompile(MOI.set, (Optimizer, MOI.ObjectiveSense, MOI.OptimizationSense))
    Base.precompile(
        MOI.set,
        (
            Optimizer,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
            MOI.ScalarAffineFunction{Float64},
        ),
    )

    functions = (MOI.ScalarAffineFunction{Float64}, MOI.SingleVariable)
    sets = (
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
    for F in functions, S in sets
        Base.precompile(MOI.add_constraint, (Optimizer, F, S))
        Base.precompile(
            MOI.get,
            (Optimizer, MOI.ConstraintPrimal, MOI.ConstraintIndex{F,S}),
        )
        Base.precompile(MOI.get, (Optimizer, MOI.ConstraintDual, MOI.ConstraintIndex{F,S}))
        Base.precompile(
            MOI.get,
            (Optimizer, MOI.ConstraintFunction, MOI.ConstraintIndex{F,S}),
        )
        Base.precompile(
            MOI.set,
            (Optimizer, MOI.ConstraintFunction, MOI.ConstraintIndex{F,S}, F),
        )
        Base.precompile(MOI.get, (Optimizer, MOI.ConstraintSet, MOI.ConstraintIndex{F,S}))
        Base.precompile(
            MOI.set,
            (Optimizer, MOI.ConstraintSet, MOI.ConstraintIndex{F,S}, S),
        )
        Base.precompile(MOI.is_valid, (Optimizer, MOI.ConstraintIndex{F,S}))
        Base.precompile(MOI.delete, (Optimizer, MOI.ConstraintIndex{F,S}))
    end
    Base.precompile(MOI.add_constraint, (Optimizer, MOI.SingleVariable, MOI.ZeroOne))
    Base.precompile(MOI.add_constraint, (Optimizer, MOI.SingleVariable, MOI.Integer))
    Base.precompile(MOI.optimize!, (Optimizer,))

    for attr in (
        MOI.TerminationStatus,
        MOI.ObjectiveValue,
        MOI.PrimalStatus,
        MOI.DualStatus,
        MOI.CallbackNodeStatus{CallbackData},
    )
        Base.precompile(MOI.get, (Optimizer, attr))
    end

    for callback in (
        CallbackFunction,
        MOI.LazyConstraintCallback,
        MOI.UserCutCallback,
        MOI.HeuristicCallback,
    )
        Base.precompile(MOI.set, (Optimizer, callback, Function))
    end

    Base.precompile(
        MOI.submit,
        (
            Optimizer,
            MOI.HeuristicSolution{CallbackData},
            Vector{MOI.VariableIndex},
            Vector{Float64},
        ),
    )
    for S in sets
        Base.precompile(
            MOI.submit,
            (Optimizer, MOI.UserCut{CallbackData}, MOI.ScalarAffineFunction{Float64}, S),
        )
        Base.precompile(
            MOI.submit,
            (
                Optimizer,
                MOI.LazyConstraint{CallbackData},
                MOI.ScalarAffineFunction{Float64},
                S,
            ),
        )
    end
end
