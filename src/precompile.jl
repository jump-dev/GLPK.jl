# Copyright (c) 2012 GLPK.jl contributors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the Licence, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(MathOptInterface.copy_to),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}})   # time: 0.1371813
    Base.precompile(Tuple{Type{Optimizer}})   # time: 0.08901196
    Base.precompile(MOI.set, (Optimizer, MOI.Silent, Bool))
    Base.precompile(MOI.add_variables, (Optimizer, Int))
    Base.precompile(MOI.add_variable, (Optimizer,))
    Base.precompile(MOI.delete, (Optimizer, MOI.VariableIndex))
    Base.precompile(MOI.get, (Optimizer, MOI.VariablePrimal, MOI.VariableIndex))

    Base.precompile(
        MOI.set,
        (Optimizer, MOI.ObjectiveSense, MOI.OptimizationSense),
    )
    Base.precompile(
        MOI.set,
        (
            Optimizer,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
            MOI.ScalarAffineFunction{Float64},
        ),
    )

    functions = (MOI.ScalarAffineFunction{Float64}, MOI.VariableIndex)
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
        Base.precompile(
            MOI.get,
            (Optimizer, MOI.ConstraintDual, MOI.ConstraintIndex{F,S}),
        )
        Base.precompile(
            MOI.get,
            (Optimizer, MOI.ConstraintFunction, MOI.ConstraintIndex{F,S}),
        )
        Base.precompile(
            MOI.set,
            (Optimizer, MOI.ConstraintFunction, MOI.ConstraintIndex{F,S}, F),
        )
        Base.precompile(
            MOI.get,
            (Optimizer, MOI.ConstraintSet, MOI.ConstraintIndex{F,S}),
        )
        Base.precompile(
            MOI.set,
            (Optimizer, MOI.ConstraintSet, MOI.ConstraintIndex{F,S}, S),
        )
        Base.precompile(MOI.is_valid, (Optimizer, MOI.ConstraintIndex{F,S}))
        Base.precompile(MOI.delete, (Optimizer, MOI.ConstraintIndex{F,S}))
    end
    Base.precompile(
        MOI.add_constraint,
        (Optimizer, MOI.VariableIndex, MOI.ZeroOne),
    )
    Base.precompile(
        MOI.add_constraint,
        (Optimizer, MOI.VariableIndex, MOI.Integer),
    )
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
            (
                Optimizer,
                MOI.UserCut{CallbackData},
                MOI.ScalarAffineFunction{Float64},
                S,
            ),
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
