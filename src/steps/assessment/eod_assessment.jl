using .Networks

using Dates
using JuMP
using Printf
using Parameters

@with_kw mutable struct EODAssessmentConfigs <: AbstractRunnableConfigs
    cut_prod_coeff = 100 #must be greater than LoL coeff in objective
    inj_prod_coeff = 10 #must be greater than LoL coeff in objective
    loss_of_load_coeff = 1 #LoL
    out_path = nothing
    problem_name = "EODAssessment"
end

@with_kw struct EODAssessment <: AbstractRunnable
    configs::EODAssessmentConfigs = EODAssessmentConfigs()
end

function run(runnable::EODAssessment, ech::Dates.DateTime, TS::Vector{Dates.DateTime}, context::AbstractContext)
    problem_name_l = @sprintf("eod_assessment_%s", ech)
    runnable.configs.out_path = context.out_dir
    runnable.configs.problem_name = problem_name_l

    # @assert( runnable.configs.loss_of_load_penalty >=
    #         biggest_min_imposition(get_tso_actions(context)) * runnable.configs.cut_prod_penalty)
    @assert( runnable.configs.loss_of_load_coeff
                <= runnable.configs.inj_prod_coeff
                    <= runnable.configs.cut_prod_coeff )

    #FIXME : only create vars for started units

    return eod_assessment(get_network(context),
                        TS,
                        get_assessment_uncertainties(context),
                        get_tso_actions(context),
                        runnable.configs
                    )
end

function eod_assessment(network,
                        TS,
                        assessment_uncertainties,
                        tso_actions,
                        configs::EODAssessmentConfigs
                    )
    model_container_l = formulate_eod_assessment(network, TS, assessment_uncertainties, tso_actions, configs)
    solve!(model_container_l, configs.problem_name, configs.out_path)

    for var in all_variables(model_container_l.model)
        println(name(var), " = ", value(var))
    end

    return model_container_l
end

function biggest_min_imposition(tso_actions::TSOActions)
    return maximum(x->x[2][1], get_impositions(tso_actions))
end

