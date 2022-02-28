try
    using Xpress;
    global OPTIMIZER = Xpress.Optimizer
catch e_xpress
    if isa(e_xpress, ArgumentError)
        try
            using CPLEX;
            global OPTIMIZER = CPLEX.Optimizer
        catch e_cplex
            if isa(e_cplex, ArgumentError)
                using Cbc;
                global OPTIMIZER = Cbc.Optimizer
            else
                throw(e_cplex)
            end
        end
    else
        throw(e_xpress)
    end
end

"""
Possible status values for a pscopf model container

    - pscopf_OPTIMAL : a solution that does not use slacks was retrieved
    - pscopf_INFEASIBLE : no solution was retrieved
    - pscopf_FEASIBLE : non-optimal solution was retrieved
    - pscopf_UNSOLVED : model is not solved yet
"""
@enum PSCOPFStatus begin
    pscopf_OPTIMAL
    pscopf_INFEASIBLE
    pscopf_FEASIBLE
    pscopf_UNSOLVED
end


abstract type AbstractModelContainer end

function get_model(model_container::AbstractModelContainer)::Model
    return model_container.model
end

function get_status(model_container_p::AbstractModelContainer)::PSCOPFStatus
    solver_status_l = termination_status(get_model(model_container_p))

    if solver_status_l == OPTIMIZE_NOT_CALLED
        return pscopf_UNSOLVED
    elseif solver_status_l == INFEASIBLE
        @error "model status is infeasible!"
        return pscopf_INFEASIBLE
    elseif solver_status_l == OPTIMAL
        return pscopf_OPTIMAL
    else
        @warn "solver termination status was not optimal : $(solver_status_l)"
        return pscopf_FEASIBLE
    end
end

function solve!(model_container::AbstractModelContainer,
                problem_name="problem", out_folder=".",
                optimizer=OPTIMIZER)
    mkpath(out_folder)

    model_l = get_model(model_container)
    set_optimizer(model_l, optimizer);

    model_file_l = joinpath(out_folder, problem_name*".lp")
    write_to_file(model_l, model_file_l)

    log_file_l = joinpath(out_folder, problem_name*".log")
    redirect_to_file(log_file_l) do
        optimize!(model_l)
    end
end






abstract type AbstractGeneratorModel end
abstract type AbstractImposableModel <: AbstractGeneratorModel end
abstract type AbstractLimitableModel <: AbstractGeneratorModel end

abstract type AbstractObjectiveModel end


###################################
# Utils
###################################

"""
    redirect_to_file(f::Function, file_p::String, mode_p="w")

Execute function `f` while redirecting C and Julia level stdout to the file file_p.
Note that `file_p` is open in write mode by default.

# Arguments
- `f::Function` : the function to execute
- `file_p::String` : name of the file to print to
- `mode_p` : open mode of the file (defaults to "w")
"""
function redirect_to_file(f::Function, file_p::String, mode_p="w")
    open(file_p, mode_p) do file_l
        Base.Libc.flush_cstdio()
        redirect_stdout(file_l) do
            f()
            Base.Libc.flush_cstdio()
        end
    end
end
