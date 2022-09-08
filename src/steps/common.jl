using .Networks

using JuMP
using Dates
using MathProgBase

try
    using Xpress;
    global OPTIMIZER = Xpress.Optimizer
    global OPTIMIZER_NAME = "Xpress"
catch e_xpress
    if isa(e_xpress, ArgumentError)
        try
            using CPLEX;
            global OPTIMIZER = CPLEX.Optimizer
            global OPTIMIZER_NAME = "CPLEX"
        catch e_cplex
            if isa(e_cplex, ArgumentError)
                using Cbc;
                global OPTIMIZER = Cbc.Optimizer
                global OPTIMIZER_NAME = "Cbc"
            else
                throw(e_cplex)
            end
        end
    else
        throw(e_xpress)
    end
end
@info("used optimizer: ", OPTIMIZER)

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
    pscopf_HAS_SLACK
    pscopf_UNSOLVED
end


abstract type AbstractModelContainer end

abstract type AbstractGeneratorModel end
abstract type AbstractPilotableModel <: AbstractGeneratorModel end
abstract type AbstractLimitableModel <: AbstractGeneratorModel end

abstract type AbstractLoLModel end

abstract type AbstractObjectiveModel end

# TODO : use a proper struct for by scenario variables
# eg: UncertainValue{VariableRef}
# or maybe a similar other struct cause depending on "ech"
#    we either will need by scenario vars (eg. injection/commitment before DP/DMO for pilotables)
#    or we will need one variable (eg. injection at or after DP for pilotables)
# if need a firm value (at/after DP or DMO) call a link_scenarios(::AbstractModel, ::)
# which adds @constraint(model, by_scenario_vars[s] == firm_variable)
# or maybe no need to if we create a single variable for scenarios right from the beginning
# and have proper getters too get_var(::, s) -> scenario's var or missing
# and have proper getters too get_var(::) -> firm var or error

struct BilevelModelContainer{U,L,K} <: AbstractModelContainer
    model::Model
    upper::U
    lower::L
    kkt_model::K
end

# AbstractModelContainer
###########################

function get_model(model_container::AbstractModelContainer)::AbstractModel
    return model_container.model
end

function get_status(model_container_p::AbstractModelContainer)::PSCOPFStatus
    solver_status_l = termination_status(get_model(model_container_p))

    if solver_status_l == MOI.OPTIMIZE_NOT_CALLED
        return pscopf_UNSOLVED
    elseif solver_status_l == MOI.INFEASIBLE
        return pscopf_INFEASIBLE
    elseif solver_status_l == MOI.OPTIMAL
        if has_positive_slack(model_container_p)
            return pscopf_HAS_SLACK
        else
            return pscopf_OPTIMAL
        end
    else
        return pscopf_FEASIBLE
    end
end


"""
    set_start_values!(model)
Set start values for the model to the solved values in the input model.
Note :
    The model is no longer optimized, it goes back to a OPTIMIZE_NOT_CALLED state
"""
function set_start_values!(model::AbstractModel)
    vars = all_variables(model)
    sol_vals = value.(vars)
    set_start_value.(vars, sol_vals)

    return model
end

function get_objective_model(model_container::AbstractModelContainer)
    return model_container.objective_model
end

function solve_step1!(model_container::AbstractModelContainer, configs::AbstractRunnableConfigs)
    model_l = get_model(model_container)

    # delete old bound if needed
    deltas_cstr = get_deltas_bounding_constraint(model_container)
    if !ismissing(deltas_cstr)
        delete(model_l, deltas_cstr)
        unregister(model_l, :deltas_cstr)
        set_deltas_bounding_constraint!(model_container, missing)
    end

    obj = get_objective_model(model_container).full_obj_1
    @objective(model_l, Min, obj)
    solve!(model_container, configs.problem_name*"_step1", configs.out_path)
    @info "deltas current value : $(value(get_objective_model(model_container).deltas))"
    @info "step2 objective current value : $(value(get_objective_model(model_container).full_obj_2))"
end

function solve_step2!(model_container::AbstractModelContainer, configs::AbstractRunnableConfigs)
    if (get_status(model_container)!=pscopf_INFEASIBLE
        && value(get_objective_model(model_container).deltas)>0 )

        model_l = get_model(model_container)

        # before invalidating the model,
        value_sum_deltas_l = value(get_objective_model(model_container).deltas)
        set_start_values!(model_l)

        bound_sum_p_deltas!(model_container, value_sum_deltas_l)

        obj = get_objective_model(model_container).full_obj_2
        @objective(model_l, Min, obj)
        solve!(model_container, configs.problem_name*"_step2", configs.out_path)
        @info "step 1 objective current value (deltas+penalty) : $(value(get_objective_model(model_container).full_obj_1))"
        @info "step 1 objective current value (deltas) : $(value(get_objective_model(model_container).deltas))"
    end
end

function solve_2steps_deltas!(model_container::AbstractModelContainer, configs::AbstractRunnableConfigs)
    solve_step1!(model_container, configs)
    solve_step2!(model_container, configs)
end

function bound_sum_p_deltas!(model_container::AbstractModelContainer, value_sum_deltas::Float64)
    model_l = get_model(model_container)
    deltas_expr = get_objective_model(model_container).deltas

    cstr = @constraint(model_l, deltas_expr<=value_sum_deltas)
    set_deltas_bounding_constraint!(model_container, cstr)
    return model_container
end


function solve!(model::AbstractModel,
                problem_name="problem", out_folder=nothing,
                optimizer=OPTIMIZER)
    problem_name_l = valid_filename(problem_name)
    set_optimizer(model, optimizer);

    if !isnothing(out_folder)
        mkpath(out_folder)
        if get_config("LP_FILES")
            @timeit TIMER_TRACKS "write_lp" begin

                model_file_l = joinpath(out_folder, problem_name_l*".lp")
                write_to_file(model, model_file_l)

                model_file_l_new = joinpath(out_folder, problem_name_l*"_new.lp")

                MOI.Utilities.attach_optimizer(model)
                internalModel = unsafe_backend(model)
                if OPTIMIZER_NAME == "Xpress"
                    Xpress.writeprob(internalModel.inner, model_file_l_new, "-l");
                elseif OPTIMIZER_NAME == "CPLEX"
                    CPLEX.CPXwriteprob(internalModel.env, internalModel.lp, model_file_l_new, "lp")
                end

                @info "LP files written."
            end
        end

        log_file_l = joinpath(out_folder, problem_name_l*".log")
    else
        log_file_l = devnull
    end

    time_limit_l = get_config("PSCOPF_TIME_LIMIT_IN_SECONDS")
    if !isnothing(time_limit_l)
        set_time_limit_sec(model, time_limit_l)
    end
    if get_config("PSCOPF_REDIRECT_LOG")
        redirect_to_file(log_file_l) do
            optimize!(model)
        end
    else
        optimize!(model)
    end

    return model
end

function solve!(model_container::AbstractModelContainer, problem_name, out_path)
    model_l = get_model(model_container)

    if get_config("FIX_NULL_LOL")
        freeze_null_lol!(model_container)
    end

    @info @sprintf("solve!(%s)", problem_name)
    solve!(model_l, problem_name, out_path)

    solve_pscopf_status = get_status(model_container)
    @info "pscopf model solve time: $(solve_time(model_l))"
    @info "pscopf model status: $(solve_pscopf_status)"
    @info "Termination status : $(termination_status(model_l))"
    @info "Objective value : $(objective_value(model_l))"
    if solve_pscopf_status == pscopf_UNSOLVED
        @error "Model not solved!"
    elseif solve_pscopf_status == pscopf_INFEASIBLE
        @error "model status is infeasible!"
    elseif solve_pscopf_status == pscopf_HAS_SLACK
        @warn "model solved optimally but slack variables were used!"
    end

    total_lol_l = total_lol(model_container)
    if !isnothing(total_lol_l)
        @info @sprintf("LoL = %.10f", total_lol_l)
    end

    return model_container
end

function write_solution(filename, model)
    open(filename, "w") do file_l
        for var in all_variables(model)
            Base.write(file_l, @sprintf("%60s = %f\n", name(var), value(var)))
        end
    end
end

function has_positive_slack(model_container)::Bool
    error("unimplemented")
end

function has_global_lol(model_container::AbstractModelContainer)
    return hasproperty(model_container.lol_model, :p_global_loss_of_load)
end

function requires_linking(firmness::DecisionFirmness, do_link::Bool=false)::Bool
    return do_link || (firmness in [DECIDED, TO_DECIDE])
end


function get_deltas_bounding_constraint(model_container::AbstractModelContainer)::Union{ConstraintRef, Missing}
    return model_container.deltas_bounding_constraint
end

function set_deltas_bounding_constraint!(model_container::AbstractModelContainer, cstr::Union{ConstraintRef,Missing})
    model_container.deltas_bounding_constraint = cstr
end

# AbstractGeneratorModel
############################

function sum_injections(generator_model::AbstractGeneratorModel,
                        ts::Dates.DateTime, s::String)
    sum_l = 0.
    for ((_,ts_l,s_l), var_l) in get_p_injected(generator_model)
        if (ts_l,s_l) == (ts, s)
            sum_l += var_l
        end
    end
    return sum_l
end

function get_p_injected(model_container::AbstractModelContainer, type::Networks.GeneratorType)
    if type == Networks.LIMITABLE
        return model_container.limitable_model.p_injected
    elseif type == Networks.PILOTABLE
        return model_container.pilotable_model.p_injected
    end
    return nothing
end

function get_p_injected(generator_model::AbstractGeneratorModel)
    return generator_model.p_injected
end

function has_injections(generator_model::AbstractGeneratorModel)
    return hasproperty(generator_model, :p_injected)
end

function add_injection_vars!(model::AbstractModel, generator_model::AbstractGeneratorModel,
                            generators_list::Vector{Networks.Generator}, target_timepoints, scenarios,
                            prefix::String)
    for gen in generators_list
        gen_id = Networks.get_id(gen)
        for ts in target_timepoints
            for s in scenarios
                add_p_injected!(generator_model, model,
                                gen_id, ts, s,
                                prefix)
            end
        end
    end
end

function add_p_injected!(generator_model::AbstractGeneratorModel, model::AbstractModel,
                        gen_id::String, ts::DateTime, s::String,
                        prefix::String
                        )::AbstractVariableRef
    name =  @sprintf("%sP_injected[%s,%s,%s]", prefix, gen_id, ts, s)
    var_l = get_p_injected(generator_model)[gen_id, ts, s] = @variable(model, base_name=name, lower_bound=0.)
    return var_l
end

function add_p_delta!(generator_model::AbstractGeneratorModel, model::Model,
                        gen_id::String, ts::DateTime, s::String,
                        p_reference::Float64,
                        prefix::String)::VariableRef
    @assert(has_injections(generator_model))
    p_injected = get_p_injected(generator_model)
    deltas = generator_model.delta_p

    name =  @sprintf("%sDelta_p[%s,%s,%s]", gen_id, ts, s, prefix)
    deltas[gen_id, ts, s] = @variable(model, base_name=name, lower_bound=0.)
    @constraint(model, deltas[gen_id, ts, s] >= p_injected[gen_id, ts, s] - p_reference)
    @constraint(model, deltas[gen_id, ts, s] >= p_reference - p_injected[gen_id, ts, s])

    return generator_model.delta_p[gen_id, ts, s]
end

function add_delta_p_vars!(model::AbstractModel, generator_model::AbstractGeneratorModel,
                        generators_list::Vector{Networks.Generator}, target_timepoints, scenarios,
                        reference_market_schedule::Schedule,
                        prefix::String)
    for gen in generators_list
        gen_id = Networks.get_id(gen)
        for ts in target_timepoints
            for s in scenarios
                p_ref = get_prod_value(reference_market_schedule, gen_id, ts, s)
                p_ref = ismissing(p_ref) ? 0. : p_ref
                add_p_delta!(generator_model, model, gen_id, ts, s, p_ref, prefix)
            end
        end
    end
end


# AbstractLimitableModel
############################

function add_limitables_vars!(model_container::AbstractModelContainer, target_timepoints, scenarios,
                            limitables_list::Union{Missing,Vector{Networks.Generator}}=missing,
                            reference_market_schedule::Union{Missing,Schedule}=missing;
                            global_capping_vars::Bool=false, global_injection_vars::Bool=false,
                            local_capping_vars::Bool=false, injection_vars::Bool=false, limit_vars::Bool=false, delta_vars::Bool=false,
                            prefix::String=""
                            )
    model = get_model(model_container)
    limitable_model = model_container.limitable_model

    if injection_vars
        add_injection_vars!(model, limitable_model,
                          limitables_list, target_timepoints, scenarios, prefix)

        if delta_vars
            if ismissing(reference_market_schedule)
                error("reference_market_schedule argument is mandatory to define delta variables!")
            end
            add_delta_p_vars!(model, limitable_model,
                            limitables_list, target_timepoints, scenarios, reference_market_schedule, prefix)
        end
    end

    if limit_vars
        add_limitation_vars!(model, limitable_model, limitables_list, target_timepoints, scenarios, prefix)
    end

    if global_capping_vars
        add_global_capping_vars!(model, limitable_model, target_timepoints, scenarios, prefix)
    end

    if global_injection_vars
        add_global_injection_vars!(model, limitable_model, target_timepoints, scenarios, prefix)
    end

    if local_capping_vars
        add_local_capping_vars!(model, limitable_model, limitables_list, target_timepoints, scenarios, prefix)
    end

end


function sum_capping(limitable_model::AbstractLimitableModel, ts,s)
    return get_global_capping(limitable_model)[ts,s]
end
function sum_capping(limitable_model::AbstractLimitableModel, ts, s, ::Networks.Network)
    return get_global_capping(limitable_model)[ts,s]
end


function get_global_injection(limitable_model::AbstractLimitableModel)
    return limitable_model.p_global_injected
end
function get_global_capping(limitable_model::AbstractLimitableModel)
    return limitable_model.p_global_capping
end
function get_local_capping(limitable_model::AbstractLimitableModel)
    return limitable_model.p_capping
end


function add_global_capping_vars!(model::AbstractModel, limitable_model::AbstractLimitableModel, target_timepoints, scenarios,
                                prefix::String)
    for ts in target_timepoints
        for s in scenarios
            name =  @sprintf("%sP_global_capping[%s,%s]", prefix, ts, s)
            get_global_capping(limitable_model)[ts, s] = @variable(model, base_name=name, lower_bound=0.)
        end
    end
end


function add_global_injection_vars!(model::AbstractModel, limitable_model::AbstractLimitableModel, target_timepoints, scenarios,
                                prefix::String)
    for ts in target_timepoints
        for s in scenarios
            name =  @sprintf("%sP_global_injected[%s,%s]", prefix, ts, s)
            get_global_injection(limitable_model)[ts, s] = @variable(model, base_name=name, lower_bound=0.)
        end
    end
end


function add_local_capping_vars!(model::AbstractModel, limitable_model::AbstractLimitableModel, limitables_list, target_timepoints, scenarios,
                                prefix::String)
    for gen in limitables_list
        gen_id = Networks.get_id(gen)
        for ts in target_timepoints
            for s in scenarios
                name =  @sprintf("%sP_capping[%s,%s,%s]", prefix, gen_id, ts, s)
                get_local_capping(limitable_model)[gen_id, ts, s] = @variable(model, base_name=name, lower_bound=0.)
            end
        end
    end
end


function global_capping_constraints!(model::AbstractModel,
                                    limitable_model::AbstractLimitableModel, limitables_list,
                                    target_timepoints, scenarios,
                                    uncertainties_at_ech::UncertaintiesAtEch;
                                    min_cap::SortedDict{Tuple{DateTime,String},V}=SortedDict{Tuple{DateTime,String},Float64}(),
                                    tso_actions::TSOActions=TSOActions(),
                                    prefix::String="") where V <: Union{AbstractVariableRef,Float64}
    global_capping_vars = get_global_capping(limitable_model)
    for ts in target_timepoints
        for s in scenarios
            prod_capacity = compute_prod(uncertainties_at_ech, limitables_list, ts, s)
            c_name = @sprintf("%sc_max_global_capped[%s,%s]",prefix,ts,s)
            @constraint(model, global_capping_vars[ts, s] <= prod_capacity , base_name = c_name)

            limitations_capped = compute_capped(uncertainties_at_ech, get_limitations(tso_actions), limitables_list, ts, s)
            c_name = @sprintf("%sc_global_capped_by_limitations[%s,%s]",prefix,ts,s)
            @constraint(model, global_capping_vars[ts, s] >= limitations_capped , base_name = c_name)

            if haskey(min_cap, (ts,s))
                c_name = @sprintf("%sc_min_global_capping[%s,%s]",prefix,ts,s)
                @constraint(model, global_capping_vars[ts, s] >= min_cap[ts,s] , base_name = c_name)
            end
        end
    end
end

function global_injected_constraints!(model::AbstractModel,
                                    limitable_model::AbstractLimitableModel, limitables_list,
                                    target_timepoints, scenarios,
                                    uncertainties_at_ech::UncertaintiesAtEch;
                                    prefix::String="")
    global_injected_vars = get_global_injection(limitable_model)
    for ts in target_timepoints
        for s in scenarios
            prod_capacity = compute_prod(uncertainties_at_ech, limitables_list, ts, s)
            c_name = @sprintf("%sc_max_global_limitables_injection[%s,%s]",prefix,ts,s)
            @constraint(model, global_injected_vars[ts, s] <= prod_capacity , base_name = c_name)
        end
    end
end

function add_limitation_vars!(model::AbstractModel,
                            limitable_model::AbstractLimitableModel,
                            limitables_list::Vector{Networks.Generator}, target_timepoints, scenarios,
                            prefix::String)
    p_limit = limitable_model.p_limit
    b_is_limited = limitable_model.b_is_limited
    p_limit_x_is_limited = limitable_model.p_limit_x_is_limited

    for limitable_gen in limitables_list
        gen_id = Networks.get_id(limitable_gen)
        pmax = Networks.get_p_max(limitable_gen)
        for ts in target_timepoints
            for s in scenarios
                name =  @sprintf("%sP_limit[%s,%s,%s]", prefix, gen_id, ts, s)
                p_limit[gen_id, ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=pmax)

                name =  @sprintf("%sB_is_limited[%s,%s,%s]", prefix, gen_id, ts, s)
                b_is_limited[gen_id, ts, s] = @variable(model, base_name=name, binary=true)

                name =  @sprintf("%sP_limit_x_is_limited[%s,%s,%s]", prefix, gen_id, ts, s)
                p_limit_x_is_limited[gen_id, ts, s] = add_prod_vars!(model,
                                                                    p_limit[gen_id, ts, s],
                                                                    b_is_limited[gen_id, ts, s],
                                                                    pmax,
                                                                    name
                                                                    )
            end
        end
    end
end

function bound_limit!(model::AbstractModel,
                    limitable_model::AbstractLimitableModel, limitables_list_l, target_timepoints, scenarios)
    for gen in limitables_list_l
        gen_id = Networks.get_id(gen)
        for ts in target_timepoints
            for s in scenarios
                p_enr = Networks.get_p_max(gen)
                c_name = @sprintf("c_limitable_max[%s,%s,%s]",gen_id,ts,s)
                @constraint(model, limitable_model.p_limit[gen_id,ts,s] <= p_enr, base_name=c_name)
            end
        end
    end
end
function inject_at_limit!(model::AbstractModel,
                        limitable_model::AbstractLimitableModel,
                        limitables_list_l, target_timepoints, scenarios,
                        uncertainties_at_ech::UncertaintiesAtEch)
    p_limit = limitable_model.p_limit
    b_is_limited = limitable_model.b_is_limited
    p_limit_x_is_limited = limitable_model.p_limit_x_is_limited
    for gen in limitables_list_l
        gen_id = Networks.get_id(gen)
        for ts in target_timepoints
            for s in scenarios
                injection_var = get_p_injected(limitable_model)[gen_id,ts,s]
                p_enr = min(Networks.get_p_max(gen), get_uncertainties(uncertainties_at_ech, gen_id, ts, s))
                #inj[g,ts,s] = min{p_limit[g,ts,s], uncertainties(g,ts,s), pmax(g)}
                name = @sprintf("pinj_lim_by_uncertainty[%s,%s,%s]",gen_id,ts,s)
                @constraint(model, injection_var <= p_enr, base_name = name)
                name = @sprintf("pinj_lim_by_limit[%s,%s,%s]",gen_id,ts,s)
                @constraint(model, injection_var <= p_limit[gen_id,ts,s], base_name = name)
                name = @sprintf("pinj_at_limit_or_uncertainty[%s,%s,%s]",gen_id,ts,s)
                @constraint(model, injection_var ==
                                (1-b_is_limited[gen_id, ts, s]) * p_enr + p_limit_x_is_limited[gen_id, ts, s],
                            base_name = name)
            end
        end
    end
end
function limitable_power_constraints!(model::AbstractModel,
                                    limitable_model::AbstractLimitableModel,
                                    limitables_list_l, target_timepoints, scenarios,
                                    firmness, uncertainties_at_ech::UncertaintiesAtEch;
                                    always_link_scenarios::Bool=false,
                                    prefix::String="")
    bound_limit!(model,
                limitable_model, limitables_list_l, target_timepoints, scenarios)
    inject_at_limit!(model,
                    limitable_model, limitables_list_l, target_timepoints, scenarios, uncertainties_at_ech)

    for gen in limitables_list_l
        @assert(Networks.get_type(gen) == Networks.LIMITABLE)
        gen_id = Networks.get_id(gen)
        add_scenarios_linking_constraints!(model, gen,
                                        limitable_model.p_limit,
                                        target_timepoints, scenarios,
                                        get_power_level_firmness(firmness, gen_id),
                                        always_link_scenarios,
                                        prefix="limitable_power_c_linking"*prefix
                                        )
    end

end

function add_p_limit!(limitable_model::AbstractLimitableModel, model::AbstractModel,
                        gen_id::String, ts::Dates.DateTime,
                        scenarios::Vector{String},
                        pmax,
                        inject_uncertainties::InjectionUncertainties,
                        decision_firmness::DecisionFirmness, #by ts
                        always_link_scenarios=false
                        )
    b_is_limited = limitable_model.b_is_limited
    p_limit_x_is_limited = limitable_model.p_limit_x_is_limited
    p_limit = limitable_model.p_limit

    for s in scenarios
        name =  @sprintf("P_limit[%s,%s,%s]", gen_id, ts, s)
        p_limit[gen_id, ts, s] = @variable(model, base_name=name, lower_bound=0., upper_bound=pmax)

        injection_var = limitable_model.p_injected[gen_id, ts, s]

        name =  @sprintf("B_is_limited[%s,%s,%s]", gen_id, ts, s)
        b_is_limited[gen_id, ts, s] = @variable(model, base_name=name, binary=true)

        name =  @sprintf("P_limit_x_is_limited[%s,%s,%s]", gen_id, ts, s)
        p_limit_x_is_limited[gen_id, ts, s] = add_prod_vars!(model,
                                                            p_limit[gen_id, ts, s],
                                                            b_is_limited[gen_id, ts, s],
                                                            pmax,
                                                            name
                                                            )

        #inj[g,ts,s] = min{p_limit[g,ts,s], uncertainties(g,ts,s), pmax(g)}
        name = @sprintf("c1_pinj_lim[%s,%s,%s]",gen_id,ts,s)
        @constraint(model, injection_var <= p_limit[gen_id, ts, s], base_name = name)
        name = @sprintf("c2_pinj_lim[%s,%s,%s]",gen_id,ts,s)
        p_enr = min(get_uncertainties(inject_uncertainties, ts, s), pmax)
        @constraint(model, injection_var ==
                        (1-b_is_limited[gen_id, ts, s]) * p_enr + p_limit_x_is_limited[gen_id, ts, s],
                    base_name = name)
    end

    # NOTE : DECIDED here does not hold its meaning. FIRM is more expressive.
    #       p_limit can always be changed in the future horizons (it is not really decided)
    #DECIDED indicates that we are past the limitable's DP => need a common decision for all scenarios
    if requires_linking(decision_firmness, always_link_scenarios)
        link_scenarios!(model, p_limit, gen_id, ts, scenarios)
    end

    return limitable_model, model
end

function local_capping_constraints!(model::AbstractModel,
                            limitable_model::AbstractLimitableModel,
                            limitables_list, target_timepoints, scenarios,
                            uncertainties_at_ech;
                            prefix::String="")
    p_capping = get_local_capping(limitable_model)
    p_injected = get_p_injected(limitable_model)

    for gen in limitables_list
        gen_id = Networks.get_id(gen)
        for ts in target_timepoints
            for s in scenarios
                p_uncert = get_uncertainties(uncertainties_at_ech, gen_id, ts, s)
                c_name =  @sprintf("%sc_ub_p_capping[%s,%s,%s]", prefix, gen_id, ts, s)
                @constraint(model, p_capping[gen_id, ts, s] <= p_uncert, base_name=c_name)
                c_name =  @sprintf("%sc_link_injection_to_capping[%s,%s,%s]", prefix, gen_id, ts, s)
                @constraint(model, p_capping[gen_id, ts, s] == p_uncert - p_injected[gen_id,ts,s], base_name=c_name)
            end
        end
    end
end

# AbstractPilotableModel
############################

function get_b_on(pilotable_model::AbstractPilotableModel)
    return pilotable_model.b_on
end

function get_b_start(pilotable_model::AbstractPilotableModel)
    return pilotable_model.b_start
end

function get_p_imposition_min(pilotable_model::AbstractPilotableModel)
    return pilotable_model.p_imposition_min
end

function get_p_imposition_max(pilotable_model::AbstractPilotableModel)
    return pilotable_model.p_imposition_max
end


function add_pilotables_vars!(model_container::AbstractModelContainer,
                            pilotables_list, target_timepoints, scenarios,
                            reference_market_schedule::Union{Schedule,Missing}=missing;
                            injection_vars::Bool=false, commitment_vars::Bool=false,
                            delta_vars::Bool=false,
                            imposition_vars::Bool=false,
                            prefix::String="")
    model = get_model(model_container)
    pilotable_model = model_container.pilotable_model

    if injection_vars
        add_injection_vars!(model, pilotable_model,
                            pilotables_list, target_timepoints, scenarios, prefix)

        if delta_vars
            if ismissing(reference_market_schedule)
                error("reference_market_schedule argument is mandatory to define delta variables!")
            end
            add_delta_p_vars!(model, pilotable_model,
                                pilotables_list, target_timepoints, scenarios, reference_market_schedule, prefix)
        end
    end

    if commitment_vars
        add_unit_commitment_vars!(model, pilotable_model,
                                pilotables_list, target_timepoints, scenarios, prefix)
    end

    if imposition_vars
        add_imposition_vars(model, pilotable_model, pilotables_list, target_timepoints, scenarios, prefix)
    end
end


function add_unit_commitment_vars!(model::AbstractModel, pilotable_model::AbstractPilotableModel,
                                    pilotables_list, target_timepoints, scenarios,
                                    prefix::String)
    for gen in pilotables_list
        if Networks.needs_commitment(gen)
            gen_id = Networks.get_id(gen)
            for ts in target_timepoints
                for s in scenarios
                    add_b_on_start!(pilotable_model, model,
                                    gen_id, ts, s, prefix)
                end
            end
        end
    end
end

function add_imposition_vars(model::AbstractModel, pilotable_model::AbstractPilotableModel,
                            pilotables_list::Vector{Generator}, target_timepoints, scenarios,
                            prefix::String="")
    for generator in pilotables_list
        gen_id = Networks.get_id(generator)
        p_max = Networks.get_p_max(generator)
        for ts in target_timepoints
            for s in scenarios
                name =  @sprintf("%sP_imposition_min[%s,%s,%s]", prefix, gen_id, ts, s)
                get_p_imposition_min(pilotable_model)[gen_id, ts, s] = @variable(model, base_name=name,
                                                                    lower_bound=0., upper_bound=p_max)
                name =  @sprintf("%sP_imposition_max[%s,%s,%s]", prefix, gen_id, ts, s)
                get_p_imposition_max(pilotable_model)[gen_id, ts, s] = @variable(model, base_name=name,
                                                                    lower_bound=0., upper_bound=p_max)
            end
        end
    end
end

function add_b_on_start!(pilotable_model::AbstractPilotableModel, model::AbstractModel,
                    gen_id::String, ts::DateTime, s::String, prefix::String
                    )
    name =  @sprintf("%sB_on[%s,%s,%s]", prefix, gen_id, ts, s)
    get_b_on(pilotable_model)[gen_id, ts, s] = @variable(model, base_name=name, binary=true)
    name =  @sprintf("%sB_start[%s,%s,%s]", prefix, gen_id, ts, s)
    get_b_start(pilotable_model)[gen_id, ts, s] = @variable(model, base_name=name, binary=true)
end

function pilotable_power_constraints!(model::AbstractModel,
                                    pilotable_model::AbstractPilotableModel,
                                    pilotable_gen::Generator,
                                    target_timepoints,
                                    scenarios,
                                    gen_commitment_firmness::Union{Missing, SortedDict{Dates.DateTime, PSCOPF.DecisionFirmness}},
                                    gen_power_level_firmness::SortedDict{Dates.DateTime, PSCOPF.DecisionFirmness},
                                    gen_schedule::Union{Missing, GeneratorSchedule},
                                    ignore_sequencing_cstrs::Bool,
                                    always_link_scenarios::Bool,
                                    prefix::String
                                    )
    @assert(Networks.get_type(pilotable_gen) == Networks.PILOTABLE)

    gen_id = Networks.get_id(pilotable_gen)
    p_max = Networks.get_p_max(pilotable_gen)
    for ts in target_timepoints
        for s in scenarios
            c_name = @sprintf("c_generator_p_max[%s,%s,%s]",gen_id,ts,s)
            @constraint(model, get_p_injected(pilotable_model)[gen_id,ts,s] <= p_max, base_name=c_name)
        end
    end

    add_scenarios_linking_constraints!(model, pilotable_gen,
                                        get_p_injected(pilotable_model),
                                        target_timepoints, scenarios,
                                        gen_power_level_firmness,
                                        always_link_scenarios,
                                        prefix="pilotable_power_c_linking"*prefix
                                        )

    if !ignore_sequencing_cstrs
        if ismissing(gen_schedule)
            error("reference schedule is required for sequencing constraints when ignore_sequencing_cstrs=false.")
        end
        add_power_level_decided_constraints!(model,
                                            pilotable_gen,
                                            get_p_injected(pilotable_model),
                                            target_timepoints,
                                            scenarios,
                                            gen_commitment_firmness,
                                            gen_power_level_firmness,
                                            gen_schedule,
                                            )
    end
end
function pilotable_power_constraints!(model::AbstractModel,
                                    pilotable_model::AbstractPilotableModel,
                                    pilotables_list_l, target_timepoints, scenarios,
                                    firmness::Firmness,
                                    reference_schedule::Union{Missing,Schedule};
                                    ignore_sequencing_cstrs::Bool=false,
                                    always_link_scenarios::Bool=false,
                                    prefix::String="")
    if ismissing(reference_schedule) && !ignore_sequencing_cstrs
        error("reference schedule is required when ignore_sequencing_cstrs=false.")
    end
    for pilotable_gen in pilotables_list_l
        gen_id = Networks.get_id(pilotable_gen)
        pilotable_power_constraints!(model,
                        pilotable_model,
                        pilotable_gen,
                        target_timepoints,
                        scenarios,
                        get_commitment_firmness(firmness, gen_id), #can be missing for pilotables with pmin=0
                        get_power_level_firmness(firmness, gen_id),
                        ismissing(reference_schedule) ? missing : get_sub_schedule(reference_schedule, gen_id),
                        ignore_sequencing_cstrs,
                        always_link_scenarios,
                        prefix
                        )
    end
end


function unit_commitment_constraints!(model::AbstractModel,
                                pilotable_model::AbstractPilotableModel, pilotables_list_l,  target_timepoints, scenarios,
                                firmness, reference_schedule, generators_initial_state;
                                always_link_scenarios::Bool=false,
                                prefix::String="")
    for pilotable_gen in pilotables_list_l
        if !Networks.needs_commitment(pilotable_gen)
            continue
        end
        gen_id = Networks.get_id(pilotable_gen)

        add_commitment_constraints!(model,
                                    get_b_on(pilotable_model),
                                    get_b_start(pilotable_model),
                                    gen_id::String,
                                    target_timepoints,
                                    scenarios,
                                    get_initial_state(generators_initial_state, pilotable_gen))

        add_scenarios_linking_constraints!(model,
                                        pilotable_gen, get_b_on(pilotable_model),
                                        target_timepoints, scenarios,
                                        get_commitment_firmness(firmness, gen_id),
                                        always_link_scenarios,
                                        prefix="commitment_c_linking"*prefix
                                        )
        #linking b_on => linking b_start
        add_commitment_sequencing_constraints!(model, pilotable_gen,
                                            get_b_on(pilotable_model),
                                            get_b_start(pilotable_model),
                                            target_timepoints, scenarios,
                                            get_commitment_firmness(firmness, gen_id),
                                            get_sub_schedule(reference_schedule, gen_id)
                                            )

        if has_injections(pilotable_model)
            link_injection_to_commitment!(model,
                                        get_p_injected(pilotable_model), get_b_on(pilotable_model),
                                        gen_id, target_timepoints, scenarios,
                                        Networks.get_p_min(pilotable_gen), Networks.get_p_max(pilotable_gen))
        end
    end
end

function link_injection_to_commitment!(model::AbstractModel,
                                        p_injected_vars, b_on_vars,
                                        gen_id::String, target_timepoints, scenarios,
                                        p_min::Float64, p_max::Float64)
    for ts in target_timepoints
        for s in scenarios
            @constraint(model, p_injected_vars[gen_id, ts, s] <= p_max * b_on_vars[gen_id, ts, s]);
            @constraint(model, p_injected_vars[gen_id, ts, s] >= p_min * b_on_vars[gen_id, ts, s]);
        end
    end
end

function add_commitment_constraints!(model::AbstractModel,
                                    b_on_vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                    b_start_vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                    gen_id::String,
                                    target_timepoints::Vector{Dates.DateTime},
                                    scenarios::Vector{String},
                                    generator_initial_state::GeneratorState)
    for s in scenarios
        for (ts_index, ts) in enumerate(target_timepoints)
            preceding_on = (ts_index > 1) ? b_on_vars[gen_id, target_timepoints[ts_index-1], s] : float(generator_initial_state)
            @constraint(model, b_start_vars[gen_id, ts, s] <= b_on_vars[gen_id, ts, s])
            @constraint(model, b_start_vars[gen_id, ts, s] <= 1 - preceding_on)
            @constraint(model, b_start_vars[gen_id, ts, s] >= b_on_vars[gen_id, ts, s] - preceding_on)
        end
    end
    return model
end


function respect_impositions_constraints!(model::AbstractModel,
                                        injections_pilotable_model::AbstractPilotableModel,
                                        impositions_pilotable_model::AbstractPilotableModel,
                                        pilotables_list_l,  target_timepoints, scenarios;
                                        prefix::String="")
    p_injected_vars = get_p_injected(injections_pilotable_model)
    p_imposition_min = get_p_imposition_min(impositions_pilotable_model)
    p_imposition_max = get_p_imposition_max(impositions_pilotable_model)
    for pilotable_gen in pilotables_list_l
        gen_id = Networks.get_id(pilotable_gen)
        for ts in target_timepoints
            for s in scenarios
                imposition_bounds = (p_imposition_min[gen_id,ts,s], p_imposition_max[gen_id,ts,s])
                respect_impositions_constraints!(model, p_injected_vars, imposition_bounds, gen_id, ts, s, prefix)
            end
        end
    end
end
function respect_impositions_constraints!(model::AbstractModel,
                                        pilotable_model::AbstractPilotableModel, pilotables_list_l,  target_timepoints, scenarios,
                                        tso_actions::TSOActions;
                                        prefix::String="")
    p_injected_vars = get_p_injected(pilotable_model)
    for pilotable_gen in pilotables_list_l
        gen_id = Networks.get_id(pilotable_gen)
        for ts in target_timepoints
            for s in scenarios
                imposition_bounds = get_imposition(tso_actions, gen_id, ts, s)
                if !ismissing(imposition_bounds)
                    respect_impositions_constraints!(model, p_injected_vars, imposition_bounds, gen_id, ts, s, prefix*"tsoactions")
                end
            end
        end
    end
end
function respect_impositions_constraints!(model, p_injected_vars, imposition_bounds, gen_id, ts, s, prefix)
    c_name = @sprintf("%sc_min_imposition[%s,%s,%s]",prefix,gen_id,ts,s)
    @constraint(model, imposition_bounds[1] <= p_injected_vars[gen_id,ts,s], base_name=c_name)
    c_name = @sprintf("%sc_max_imposition[%s,%s,%s]",prefix,gen_id,ts,s)
    @constraint(model, p_injected_vars[gen_id,ts,s] <= imposition_bounds[2], base_name=c_name)
    @debug @sprintf("%simpositions constraints [%s,%s,%s] : [%s,%s]",
                    prefix,gen_id, ts, s, imposition_bounds[1], imposition_bounds[2])
end

function power_imposition_constraints!(model::AbstractModel, pilotable_model::AbstractPilotableModel,
                                        pilotables_list, target_timepoints, scenarios,
                                        firmness::Firmness,
                                        reference_schedule::Schedule;
                                        always_link_scenarios::Bool=false,
                                        prefix::String="")
    b_on_vars = get_b_on(pilotable_model)
    p_imposition_min = get_p_imposition_min(pilotable_model)
    p_imposition_max = get_p_imposition_max(pilotable_model)

    for generator in pilotables_list

        gen_id = Networks.get_id(generator)
        p_max = Networks.get_p_max(generator)
        p_min = Networks.get_p_min(generator)
        gen_power_firmness = get_power_level_firmness(firmness, gen_id)
        gen_commitment_firmness = get_commitment_firmness(firmness, gen_id)
        gen_schedule = get_sub_schedule(reference_schedule, gen_id)

        # Bounds
        for s in scenarios
            for ts in target_timepoints
                c_name = @sprintf("%sc_ub_p_imposition_min[%s,%s,%s]",prefix,gen_id,ts,s)
                @constraint(model, p_imposition_min[gen_id, ts, s] <= p_max, base_name=c_name);
                c_name = @sprintf("%sc_ub_p_imposition_max[%s,%s,%s]",prefix,gen_id,ts,s)
                @constraint(model, p_imposition_max[gen_id, ts, s] <= p_max, base_name=c_name);
                c_name = @sprintf("%sc_imposition_min_less_max[%s,%s,%s]",prefix,gen_id,ts,s)
                @constraint(model, p_imposition_min[gen_id, ts, s] <= p_imposition_max[gen_id, ts, s], base_name=c_name);
            end
        end

        # Firmness
        add_scenarios_linking_constraints!(model, generator,
                                    p_imposition_min,
                                    target_timepoints, scenarios,
                                    gen_power_firmness,
                                    always_link_scenarios,
                                    prefix="power_imposition_min_c_linking"*prefix
                                    )
        add_scenarios_linking_constraints!(model, generator,
                                    p_imposition_max,
                                    target_timepoints, scenarios,
                                    gen_power_firmness,
                                    always_link_scenarios,
                                    prefix="power_imposition_max_c_linking"*prefix
                                    )

        # Sequencing
        add_power_level_decided_constraints!(model,
                                    generator,
                                    p_imposition_min,
                                    target_timepoints, scenarios,
                                    gen_commitment_firmness, gen_power_firmness,
                                    gen_schedule,
                                    cstr_prefix_name="decided_p_min"
                                    )
        add_power_level_decided_constraints!(model,
                                    generator,
                                    p_imposition_max,
                                    target_timepoints, scenarios,
                                    gen_commitment_firmness, gen_power_firmness,
                                    gen_schedule,
                                    cstr_prefix_name="decided_p_max"
                                    )

        # commitment-imposition linking
        if Networks.needs_commitment(generator)
            for s in scenarios
                for ts in target_timepoints
                    # pmin B_on < P_tso_min < P_tso_max < pmax B_on
                    c_name = @sprintf("c_p_imposition_min_commitment[%s,%s,%s]",gen_id,ts,s)
                    @constraint(model, p_min * b_on_vars[gen_id, ts, s] <= p_imposition_min[gen_id, ts, s], base_name=c_name);
                    c_name = @sprintf("c_p_imposition_max_commitment[%s,%s,%s]",gen_id,ts,s)
                    @constraint(model, p_imposition_max[gen_id, ts, s] <= p_max * b_on_vars[gen_id, ts, s], base_name=c_name)
                end
            end
        end
    end
end


# AbstractLoLModel
############################

function get_global_lol(lol_model::AbstractLoLModel)
    return lol_model.p_global_loss_of_load
end
function get_global_lol(model_container::AbstractModelContainer)
    return get_global_lol(model_container.lol_model)
end

function get_local_lol(lol_model::AbstractLoLModel)
    return lol_model.p_loss_of_load
end

function get_local_lol(model_container::AbstractModelContainer)
    return get_local_lol(model_container.lol_model)
end

function has_positive_value(dict_vars::AbstractDict{T,V}) where T where V <: AbstractVariableRef
    return any(e -> value(e[2]) > get_config("lol_eps"), dict_vars)
end

function add_lol_vars!(model_container::AbstractModelContainer, target_timepoints, scenarios,
                    buses_list::Union{Missing,Vector{Networks.Bus}}=missing;
                    global_lol_vars::Bool=false, local_lol_vars::Bool=false,
                    prefix::String="")
    model = get_model(model_container)
    lol_model = model_container.lol_model

    if global_lol_vars
        add_global_lol_vars!(model, lol_model, target_timepoints, scenarios, prefix)
    end

    if local_lol_vars
        add_local_lol_vars!(model, lol_model, buses_list, target_timepoints, scenarios, prefix)
    end
end


function add_global_lol_vars!(model::AbstractModel, lol_model::AbstractLoLModel,
                            target_timepoints, scenarios,
                            prefix::String="")
    for ts in target_timepoints
        for s in scenarios
            name =  @sprintf("%sP_global_lol[%s,%s]", prefix, ts, s)
            get_global_lol(lol_model)[ts, s] = @variable(model, base_name=name, lower_bound=0.)
        end
    end
end

function total_lol(model_container::AbstractModelContainer)
    if has_global_lol(model_container)
        sum( value(lol_var) for (_,lol_var) in get_global_lol(model_container) )
    else
        sum( value(lol_var) for (_,lol_var) in get_local_lol(model_container) )
    end
end

function sum_lol(lol_model::AbstractLoLModel, ts, s, network=nothing)
    return get_global_lol(lol_model)[ts,s]
end

function global_lol_constraints!(model::AbstractModel,
                            lol_model::AbstractLoLModel, buses_list, target_timepoints, scenarios,
                            uncertainties_at_ech::UncertaintiesAtEch;
                            min_lol::SortedDict{Tuple{DateTime,String},V}=SortedDict{Tuple{DateTime,String},Float64}(),
                            prefix::String="") where V <: Union{AbstractVariableRef,Float64}
    global_lol_vars = get_global_lol(lol_model)
    for ts in target_timepoints
        for s in scenarios
            max_load = compute_load(uncertainties_at_ech, buses_list, ts, s)
            c_name = @sprintf("%sc_ub_global_lol[%s,%s]",prefix,ts,s)
            @constraint(model, global_lol_vars[ts, s] <= max_load , base_name = c_name)

            if haskey(min_lol, (ts,s))
                c_name = @sprintf("%sc_min_global_lol[%s,%s]",prefix,ts,s)
                @constraint(model, global_lol_vars[ts, s] >= min_lol[ts,s] , base_name = c_name)
            end
        end
    end
end


function add_local_lol_vars!(model::AbstractModel, lol_model::AbstractLoLModel,
                            buses_list::Vector{Networks.Bus}, target_timepoints, scenarios,
                            prefix::String="")
    for bus in buses_list
        bus_id = Networks.get_id(bus)
        for ts in target_timepoints
            for s in scenarios
                name =  @sprintf("%sP_local_lol[%s,%s,%s]", prefix, bus_id, ts, s)
                get_local_lol(lol_model)[bus_id, ts, s] = @variable(model, base_name=name, lower_bound=0.)
            end
        end
    end
end

function local_lol_constraints!(model::AbstractModel, lol_model::AbstractLoLModel,
                            buses_list, target_timepoints, scenarios, uncertainties_at_ech;
                            prefix::String="")
    local_lol_vars = get_local_lol(lol_model)
    for bus in buses_list
        bus_id = Networks.get_id(bus)
        for ts in target_timepoints
            for s in scenarios
                bus_load = get_uncertainties(uncertainties_at_ech, bus_id, ts, s)
                c_name = @sprintf("%sc_ub_local_lol[%s,%s,%s]",prefix,bus_id,ts,s)
                @constraint(model, local_lol_vars[bus_id, ts, s] <= bus_load , base_name = c_name)
            end
        end
    end
end

function freeze_null_lol!(model_container::AbstractModelContainer, )
    @info "fixing LoL to 0."

    lol_dict_vars = has_global_lol(model_container) ? get_global_lol(model_container) : get_local_lol(model_container)
    lol_vars = values(lol_dict_vars)
    fix.(lol_vars, 0., force=true) # deletes the variable's original bounds

    return model_container
end


# Objective
##################

function add_pilotable_start_cost!(obj_component::AffExpr,
                                b_start::AbstractDict{T,V}, network, gratis_starts) where T <: Tuple where V <: VariableRef
    for ((gen_id,ts,_), b_start_var) in b_start
        if (gen_id,ts) in gratis_starts
            @debug(@sprintf("ignore starting cost of %s at %s", gen_id, ts))
            continue
        end
        generator = Networks.get_generator(network, gen_id)
        gen_start_cost = Networks.get_start_cost(generator)
        add_to_expression!(obj_component,
                            b_start_var * gen_start_cost)
    end

    return obj_component
end

function add_prop_cost!(obj_component::AffExpr,
                                p_injected::AbstractDict{T,V}, network)  where T <: Tuple where V <: VariableRef
    for ((gen_id,_,_), p_injected_var) in p_injected
        generator = Networks.get_generator(network, gen_id)
        gen_prop_cost = Networks.get_prop_cost(generator)
        add_to_expression!(obj_component,
                            p_injected_var * gen_prop_cost)
    end

    return obj_component
end

function add_coeffxsum_cost!(obj_component::AffExpr,
                            vars_dict::AbstractDict{T,V}, coeff::Float64)  where T where V <: VariableRef
    for (_, var_l) in vars_dict
        add_to_expression!(obj_component, coeff * var_l)
    end

    return obj_component
end


# EOD Constraints
##################


function distribution_constraint!(model::AbstractModel,
                                global_vars::SortedDict{Tuple{DateTime,String},VariableRef},
                                local_vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                local_ids::Vector{String}, target_timepoints, scenarios;
                                cstr_prefix_name::String="distribute")
    if !isempty(local_ids)
        for ts in target_timepoints
            for s in scenarios
                vars_sum = sum(local_vars[id_l, ts, s]
                                for id_l in local_ids)
                c_name = @sprintf("c_%s[%s,%s]",cstr_prefix_name,ts,s)
                @constraint(model, global_vars[ts, s] == vars_sum, base_name=c_name)
            end
        end
    end
end

function eod_constraints!(model::AbstractModel, eod_constraints::SortedDict{Tuple{Dates.DateTime,String}, ConstraintRef},
                        pilotable_model::AbstractPilotableModel,
                        limitable_model::AbstractLimitableModel,
                        lol_model::AbstractLoLModel,
                        target_timepoints, scenarios,
                        uncertainties_at_ech, network;
                        prefix::String="")
    for ts in target_timepoints
        for s in scenarios

            prod = AffExpr()
            prod += sum_injections(pilotable_model, ts, s)
            if has_injections(limitable_model)
                prod += sum_injections(limitable_model, ts, s)
            else
                prod += compute_prod(uncertainties_at_ech, network, ts, s)
                prod -= sum_capping(limitable_model, ts, s, network)
            end

            load = AffExpr()
            load += compute_load(uncertainties_at_ech, network, ts, s)
            load -= sum_lol(lol_model, ts, s, network)

            c_name = @sprintf("c_eod[%s,%s]",ts,s)
            eod_constraints[ts,s] = @constraint(model, prod == load , base_name = c_name)

        end
    end
end


# RSO Constraints
##################


function get_flows(model_container::AbstractModelContainer)
    return model_container.flows
end
function get_rso_constraints(model_container::AbstractModelContainer)
    return model_container.rso_constraints
end

function add_rso_flows_exprs!(model_container::AbstractModelContainer,
                combinations_to_add, #container/Generator of Tuple{Networks.Branch,ts::DateTime,s::String,ptdfcase::String}
                uncertainties_at_ech, network::Networks.Network)
    for (branch_id, ts, s, ptdf_case) in combinations_to_add
        branch::Networks.Branch = Networks.safeget_branch(network, branch_id)
        flow_expr_l = flow_expr(model_container,
                                branch, ts, s, ptdf_case,
                                uncertainties_at_ech, network)
        add_rso_flow_expr!(model_container, flow_expr_l,
                        branch, ts, s, ptdf_case)
    end
end

function add_rso_flow_expr!(model_container::AbstractModelContainer,
                flow_expr::AffExpr,
                branch::Networks.Branch, ts::DateTime, s::String, ptdf_case::String)
    flows::SortedDict{Tuple{String,Dates.DateTime,String,String}, AffExpr} = get_flows(model_container)
    branch_id = Networks.get_id(branch)
    flows[branch_id, ts, s, ptdf_case] = flow_expr
end

function flow_expr(model_container::AbstractModelContainer,
                    branch::Networks.Branch, ts::DateTime, s::String, ptdf_case::String,
                    uncertainties_at_ech, network::Networks.Network)
    branch_id = Networks.get_id(branch)

    flow_l = get(get_flows(model_container), (branch_id, ts, s, ptdf_case), missing)
    if !ismissing(flow_l)
        return flow_l
    end

    flow_l = AffExpr()
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        ptdf = Networks.safeget_ptdf_elt(network, branch_id, bus_id, ptdf_case)

        # + injections limitables
        for gen in Networks.get_generators_of_type(bus, Networks.LIMITABLE)
            gen_id = Networks.get_id(gen)
            var_p_injected = get_p_injected(model_container, Networks.LIMITABLE)[gen_id, ts, s]
            add_to_expression!(flow_l, ptdf*var_p_injected)
        end

        # + injections pilotables
        for gen in Networks.get_generators_of_type(bus, Networks.PILOTABLE)
            gen_id = Networks.get_id(gen)
            var_p_injected = get_p_injected(model_container, Networks.PILOTABLE)[gen_id, ts, s]
            add_to_expression!(flow_l, ptdf*var_p_injected)
        end

        # - loads
        add_to_expression!(flow_l, -ptdf * get_uncertainties(uncertainties_at_ech, bus_id, ts, s))

        # + cutting loads ~ injections
        add_to_expression!(flow_l, ptdf * get_local_lol(model_container)[bus_id, ts, s])
    end

    return flow_l
end

function add_rso_constraints!(model_container::AbstractModelContainer,
                        combinations_to_add, #container of Tuple{Networks.Branch,ts::DateTime,s::String,ptdfcase::String}
                        network::Networks.Network;
                        prefix::String="")
    model = get_model(model_container)
    rso_constraints = get_rso_constraints(model_container)
    flows = get_flows(model_container)
    for (branch_id, ts, s, ptdf_case) in combinations_to_add
        branch::Networks.Branch = Networks.safeget_branch(network, branch_id)
        add_rso_constraint!(model, rso_constraints,
                        flows,
                        branch, ts, s, ptdf_case,
                        prefix=prefix)
    end
end

function add_rso_constraint!(model::AbstractModel,
                             rso_constraints::SortedDict{Tuple{String,Dates.DateTime,String,String}, Tuple{ConstraintRef,ConstraintRef}},
                            flows::SortedDict{Tuple{String,Dates.DateTime,String,String}, AffExpr},
                            branch::Networks.Branch, ts::DateTime, s::String, ptdf_case::String;
                            prefix::String="")

    branch_id = Networks.get_id(branch)
    flow_limit_l = Networks.safeget_limit(branch, ptdf_case)
    flow_expr_l = flows[branch_id, ts, s, ptdf_case]

    c_name_lb =  @sprintf("%sRSO[%s,%s,%s,%s]_lb", prefix, branch_id, ts, s, ptdf_case)
    c_name_ub =  @sprintf("%sRSO[%s,%s,%s,%s]_ub", prefix, branch_id, ts, s, ptdf_case)
    rso_constraints[branch_id, ts, s, ptdf_case] = (@constraint(model, -flow_limit_l <= flow_expr_l , base_name=c_name_lb),
                                                    @constraint(model, flow_expr_l <= flow_limit_l , base_name=c_name_ub))
end

function flow_val(model_container::AbstractModelContainer,
                branch::Networks.Branch, ts::DateTime, s::String, ptdf_case::String,
                uncertainties_at_ech, network::Networks.Network)::Float64
    branch_id = Networks.get_id(branch)

    flow_l = 0.
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        ptdf = Networks.safeget_ptdf_elt(network, branch_id, bus_id, ptdf_case)

        # + injections limitables
        for gen in Networks.get_generators_of_type(bus, Networks.LIMITABLE)
            gen_id = Networks.get_id(gen)
            val_p_injected = value(get_p_injected(model_container, Networks.LIMITABLE)[gen_id, ts, s])
            flow_l += (ptdf*val_p_injected)
        end

        # + injections pilotables
        for gen in Networks.get_generators_of_type(bus, Networks.PILOTABLE)
            gen_id = Networks.get_id(gen)
            val_p_injected = value(get_p_injected(model_container, Networks.PILOTABLE)[gen_id, ts, s])
            flow_l += (ptdf*val_p_injected)
        end

        # - loads
        flow_l -= (ptdf * get_uncertainties(uncertainties_at_ech, bus_id, ts, s))

        # + cutting loads ~ injections
        flow_l += (ptdf * value(get_local_lol(model_container)[bus_id, ts, s]))

    end

    return flow_l
end

function flow_val(branch::Networks.Branch, ts::DateTime, s::String, ptdf_case::String,
                uncertainties_at_ech, network::Networks.Network,
                p_values_lim, p_values_pil, p_values_lol, ptdf_p)::Float64
    branch_id = Networks.get_id(branch)

    flow_l = 0.
    for bus in Networks.get_buses(network)
        bus_id = Networks.get_id(bus)
        ptdf = ptdf_p[ptdf_case, branch_id, bus_id]

        # + injections limitables
        for gen in Networks.get_generators_of_type(bus, Networks.LIMITABLE)
            gen_id = Networks.get_id(gen)
            flow_l += ptdf * p_values_lim[(gen_id, ts, s)]
        end

        # + injections pilotables
        for gen in Networks.get_generators_of_type(bus, Networks.PILOTABLE)
            gen_id = Networks.get_id(gen)
            flow_l +=  ptdf * p_values_pil[(gen_id, ts, s)]
        end

        # - loads
        flow_l -= ptdf * get_uncertainties(uncertainties_at_ech, bus_id, ts, s)

        # + cutting loads ~ injections
        flow_l += ptdf * p_values_lol[(bus_id, ts, s)]

    end

    return flow_l
end

# Helpers
##################

function available_combinations(network, target_timepoints, scenarios)
    return ((branch_id,ts,s,case)
                for branch_id in keys(network.branches)
                for ts in target_timepoints
                for s in scenarios
                for case in Networks.get_cases(network) #if the network does not hold the PTDF data for a case it will be ignored
            )
end

function theoretical_nb_combinations(network, target_timepoints, scenarios)
    return length(Networks.get_branches(network)) * length(target_timepoints) * length(scenarios) * length(Networks.get_cases(network))
end

"""
# Note
    The number of active constraints is not exact since, in a MIP context,
     some constraints that have non-negative slack might be influential/active
"""
function log_flows(model_container,network,nb_potential_cstrs,out_folder,filename)
    filename_l = valid_filename(filename)
    flows::SortedDict{Tuple{String,Dates.DateTime,String,String}, AffExpr} = get_flows(model_container)
    if !isnothing(out_folder)
        log_file_l = joinpath(out_folder, filename_l*"_flows.log")
        open(log_file_l, "w") do file_l
            line_l = @sprintf("branch_id  %15s  s  ptdf_case  flow_value flow_limit\n", "ts")
            @debug line_l
            write(file_l, line_l)
            for ((branch_id, ts, s, ptdf_case), flow_expr) in sort_by_cut_branch(flows)
                flow_val = value(flow_expr)
                flow_lim = Networks.safeget_limit(get_branch(network, branch_id), ptdf_case)
                status = (abs(flow_val) >= flow_lim-1e-09) ? "ACTIVE" : "INACTIVE"
                line_l = @sprintf("%s %s %s %s %s %s %s\n",
                                branch_id, ts, s, ptdf_case, flow_val, flow_lim, status)
                write(file_l, line_l)
                @debug line_l
            end
        end
    end

    if get_config("CNT_ACTIVE_RSO_CSTRS")
        @warn "This section is time consuming deactivate it with : set_config!(\"CNT_ACTIVE_RSO_CSTRS\", false)"

        nb_active_cstrs = 0
        for (_, (cstr_lb,cstr_ub)) in get_rso_constraints(model_container)
            if !has_slack(cstr_lb) || !has_slack(cstr_ub)
                nb_active_cstrs += 1
            end
        end

        msg_l = @sprintf("at least %d/%d RSO constraints are active!", nb_active_cstrs, nb_potential_cstrs)
        @info msg_l
        open(get_config("TEMP_GLOBAL_LOGFILE"), "a") do file_l write(file_l, msg_l*"\n") end
    end
end

function sort_by_cut_branch(flows::SortedDict{Tuple{String,Dates.DateTime,String,String}, AffExpr})
    sort(collect(flows), by= x -> (x[1][4], x[1][1:3]...))
end


# Utils
##################

function link_scenarios!(model::AbstractModel, vars::AbstractDict{Tuple{String,DateTime,String},V},
                        gen_id::String, ts::DateTime, scenarios::Vector{String};
                        prefix::String="") where V<:AbstractVariableRef
    s1 = scenarios[1]
    for (s_index, s) in enumerate(scenarios)
        if s_index > 1
            cstr_l = @constraint(model, vars[gen_id, ts, s] == vars[gen_id, ts, s1]);
            if !isnothing(name)
                set_name(cstr_l, @sprintf("%s[%s,%s,%s]",prefix,gen_id,ts,s))
            end
        end
    end
    return model
end

function add_keep_off_constraint!(model::AbstractModel, b_on_vars, b_start_vars,
                                gen_id, ts, scenarios)
    for s in scenarios
        @constraint(model, b_on_vars[gen_id, ts, s] == 0)
        @constraint(model, b_start_vars[gen_id, ts, s] == 0)
    end
    return model
end

function add_commitment_sequencing_constraints!(model::AbstractModel,
                                            generator::Networks.Generator,
                                            b_on_vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                            b_start_vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                            target_timepoints::Vector{Dates.DateTime},
                                            scenarios::Vector{String},
                                            commitment_firmness::SortedDict{Dates.DateTime, DecisionFirmness}, #by ts
                                            generator_reference_schedule::GeneratorSchedule
                                            )
    gen_id = Networks.get_id(generator)
    for ts in target_timepoints
        if commitment_firmness[ts] == DECIDED
            reference_on_val = safeget_commitment_value(generator_reference_schedule, ts)
            if reference_on_val == OFF
                add_keep_off_constraint!(model, b_on_vars, b_start_vars, gen_id, ts, scenarios)
            end
        end
    end

    return model
end


function add_scenarios_linking_constraints!(model::AbstractModel,
                                                generator::Networks.Generator,
                                                vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                                target_timepoints::Vector{Dates.DateTime},
                                                scenarios::Vector{String},
                                                gen_firmness::SortedDict{Dates.DateTime, DecisionFirmness}, #by ts
                                                always_link::Bool;
                                                prefix::String=""
                                                )
    gen_id = Networks.get_id(generator)
    for ts in target_timepoints
        if requires_linking(gen_firmness[ts], always_link)
            link_scenarios!(model, vars, gen_id, ts, scenarios, prefix=prefix)
        end
    end

    return model
end

function add_power_level_decided_constraints!(model::AbstractModel,
                                                generator::Networks.Generator,
                                                p_injected_vars::SortedDict{Tuple{String,DateTime,String},VariableRef},
                                                target_timepoints::Vector{Dates.DateTime},
                                                scenarios::Vector{String},
                                                commitment_firmness::Union{Missing,SortedDict{Dates.DateTime, DecisionFirmness}}, #by ts
                                                power_level_firmness::SortedDict{Dates.DateTime, DecisionFirmness}, #by ts
                                                generator_reference_schedule::GeneratorSchedule;
                                                cstr_prefix_name::String="decide_level"
                                                )
    @assert(Networks.get_type(generator) == Networks.PILOTABLE)

    gen_id = Networks.get_id(generator)
    for ts in target_timepoints
        ref_commitment = missing
        if (!ismissing(commitment_firmness) && (commitment_firmness[ts] == DECIDED))
            ref_commitment = get_commitment_value(generator_reference_schedule, ts)
        end

        for s in scenarios
            if !ismissing(ref_commitment) && (ref_commitment == OFF)
                c_name = @sprintf("c_decided_level_off[%s,%s,%s]",gen_id,ts,s)
                @constraint(model, p_injected_vars[gen_id,ts,s] == 0., base_name=c_name)
            elseif power_level_firmness[ts]==DECIDED
                scheduled_prod = safeget_prod_value(generator_reference_schedule,ts)
                @debug @sprintf("imposed decided level[%s,%s,%s] : %s", gen_id, ts, s, scheduled_prod)
                c_name = @sprintf("c_%s[%s,%s,%s]",cstr_prefix_name,gen_id,ts,s)
                @constraint(model, p_injected_vars[gen_id,ts,s] == scheduled_prod, base_name=c_name)
            end
        end
    end

    return model
end


"""
    add_prod_vars!
adds to the model and returns a variable that represents the product expression (noted var_a_x_b):
   var_a * var_b where var_b is a binary variable, and var_a is a positive real variable bound by M
The following constraints are added to the model :
    var_a_x_b <= var_a
    var_a_x_b <= M * var_b
    M*(1 - var_b) + var_a_x_b >= var_a
"""
function add_prod_vars!(model::AbstractModel,
                        var_a::AbstractVariableRef,
                        var_binary::AbstractVariableRef,
                        M,
                        name
                        )::AbstractVariableRef
    if !is_binary(var_binary)
        throw(error("variable var_binary needs to be binary to express the product!"))
    end
    if lower_bound(var_a) < 0
        throw(error("variable var_a needs to be positive to express the product!"))
    end

    expr_var_a::AffExpr = var_a
    var_a_x_b = add_prod_expr_x_b!(model, expr_var_a, var_binary, M, name)

    return var_a_x_b
end

function add_prod_expr_x_b!(model::AbstractModel,
                        expr_a::AffExpr,
                        var_binary::AbstractVariableRef,
                        M,
                        name
                        )::AbstractVariableRef
    if !is_binary(var_binary)
        throw(error("variable var_binary needs to be binary to express the product!"))
    end
    if compute_lb(expr_a, -1) < 0
        c_name = @sprintf("c0_%s",name)
        @constraint(model, expr_a >= 0., base_name=c_name)
    end

    var_expra_x_b = @variable(model, base_name=name, lower_bound=0., upper_bound=M)
    c_name = @sprintf("c1_%s",name)
    @constraint(model, var_expra_x_b <= expr_a, base_name=c_name)
    c_name = @sprintf("c2_%s",name)
    @constraint(model, var_expra_x_b <= M * var_binary, base_name=c_name)
    c_name = @sprintf("c3_%s",name)
    @constraint(model, M*(1-var_binary) + var_expra_x_b >= expr_a, base_name=c_name)

    return var_expra_x_b
end

function add_dual_and_indicator!(kkt_model::Model, duals_dict, indicators_dict,
                                key, name,
                                is_positive::Bool)

    add_dual!(kkt_model, duals_dict, key, name, is_positive)

    indicators_dict[key] = @variable(kkt_model, binary=true, base_name="indicator_"*name)

    return duals_dict[key] , indicators_dict[key]
end

function add_dual!(kkt_model::Model, duals_dict,
                    key, name,
                    is_positive::Bool)
    if is_positive
        duals_dict[key] = @variable(kkt_model, lower_bound=0., base_name="dual_"*name)
    else
        duals_dict[key] = @variable(kkt_model, base_name="dual_"*name)
    end

    return duals_dict[key]
end

function formulate_complementarity_constraints!(model::Model,
                                                kkt_var::VariableRef, pos_cstr_expr, b_indicator, ub_kkt, ub_cstr)
    # complementarity constraints are supposed for constraints like : g(x) <= 0 => the duals should be positive
    #model should already have the constraints -pos_cstr_expr <= 0 and kkt_var >=0
    # reformulation adds the following constraints :
    # kkt_var <= M1 * b_indicator
    # -g(x) <= M2 * (1-b_indicator)
    @assert( lower_bound(kkt_var) >= 0 )

    name = "c1_complementarity_" * JuMP.name(kkt_var)
    @constraint(model, kkt_var <= ub_kkt*b_indicator, base_name=name)
    name = "c2_complementarity_" * JuMP.name(kkt_var)
    @constraint(model, pos_cstr_expr <= ub_cstr * (1-b_indicator), base_name=name)
    #This constraints should have already been added in primal-feasibility constraints : g_i(x)<=0
    name = "c3_complementarity_" * JuMP.name(kkt_var)
    @constraint(model, pos_cstr_expr >= 0, base_name=name)
end

function compute_ub(expr::AffExpr, big_m=nothing)
    expr_ub = expr.constant
    for (coeff, var) in linear_terms(expr)
        if coeff == 0
            continue
        elseif coeff > 0 && has_upper_bound(var)
            expr_ub += coeff * upper_bound(var)
        elseif coeff < 0 && has_lower_bound(var)
            expr_ub += coeff * lower_bound(var)
        else
            expr_ub = big_m
            break;
        end
    end

    if isnothing(expr_ub)
        error("need to specify bound for expression $(expr)")
    end

    return expr_ub
end

function compute_lb(expr::AffExpr, default_lb=nothing)
    expr_ub = expr.constant
    for (coeff, var) in linear_terms(expr)
        if coeff == 0
            continue
        elseif coeff > 0 && has_lower_bound(var)
            expr_ub += coeff * lower_bound(var)
        elseif coeff < 0 && has_upper_bound(var)
            expr_ub += coeff * upper_bound(var)
        else
            expr_ub = default_lb
            break;
        end
    end

    if isnothing(expr_ub)
        error("need to specify bound for expression $(expr)")
    end

    return expr_ub
end

function has_slack(cstr::ConstraintRef, eps=1e-09)
    has_slack_l = (normalized_rhs(cstr) - value(cstr)) >= eps
    return has_slack_l
end

