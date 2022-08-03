using .Networks

using Dates
using JuMP
using Printf
using Parameters


@with_kw mutable struct RSOAssessmentConfigs <: AbstractRunnableConfigs
    BIG_M = get_config("big_m_value") #max supposed overflow

    out_path = nothing
    problem_name = "RSOAssessment"
end

@with_kw struct RSOAssessment <: AbstractRunnable
    configs::RSOAssessmentConfigs = RSOAssessmentConfigs()
end

function run(runnable::RSOAssessment, ech::Dates.DateTime, TS::Vector{Dates.DateTime}, context::AbstractContext)
    problem_name_l = @sprintf("rso_assessment_%s", ech)
    runnable.configs.out_path = context.out_dir
    runnable.configs.problem_name = problem_name_l

    #FIXME : only create vars for started units

    return rso_assessment(get_network(context),
                        TS,
                        get_assessment_uncertainties(context),
                        get_tso_actions(context),
                        runnable.configs
                    )
end

function rso_assessment(network,
                        TS,
                        assessment_uncertainties,
                        tso_actions,
                        configs::RSOAssessmentConfigs
                    )
    model_container_l = formulate_rso_assessment(network, TS, assessment_uncertainties, tso_actions, configs)
    solve!(model_container_l, configs.problem_name, configs.out_path)

    for var in all_variables(model_container_l.model)
        println(name(var), " = ", value(var))
    end

    return model_container_l
end

