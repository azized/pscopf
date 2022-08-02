using .Networks

using Dates

################################################################################
####       TSO
################################################################################

"""
utilisé pour le mode 3:
Prend des décisions fermes vu que c'est la dernière execution du TSO
Décide de la reserve
"""
struct TSOAtFOBiLevel <: AbstractTSO
end
function run(runnable::TSOAtFOBiLevel, ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime}, context::AbstractContext)
    error("unimplemented")
    return #result
end

#### TSO COMMON :
function update_tso_schedule!(context::AbstractContext, ech, result, firmness,
                            runnable::AbstractTSO)
    tso_schedule = get_tso_schedule(context)
    tso_schedule.decider_type = DeciderType(runnable)
    tso_schedule.decision_time = ech
    error("unimplemented")
end
function update_tso_actions!(context::AbstractContext, ech, result, firmness,
                            runnable::AbstractTSO)
    error("unimplemented")
end

################################################################################
####       Firmness
################################################################################

function compute_firmness(runnable::AbstractRunnable,
                    ech::Dates.DateTime, next_ech::Union{Nothing,Dates.DateTime},
                    TS::Vector{Dates.DateTime}, context::AbstractContext)
    @debug "compute decisions firmness"
    generators = collect(Networks.get_generators(get_network(context)))
    return compute_firmness(ech, next_ech, TS, generators)
end

"""
    All decisions are Firm (DECIDED or TO_DECIDE)
"""
function compute_firmness(runnable::Union{EnergyMarketAtFO,TSOAtFOBiLevel},
                    ech::Dates.DateTime, next_ech::Union{Nothing,Dates.DateTime},
                    TS::Vector{Dates.DateTime}, context::AbstractContext)
    @debug "compute decisions firmness : only firm decisions (DECIDED or TO_DECIDE)"
    next_ech = nothing
    generators = collect(Networks.get_generators(get_network(context)))
    return compute_firmness(ech, next_ech, TS, generators)
end


################################################################################
####       Utils
################################################################################

struct Assessment <: AbstractAssessment
end
function run(runnable::Assessment, ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime}, context::AbstractContext)
    return #result
end

struct EnterFO <: AbstractRunnable
end
function run(runnable::EnterFO, ech::Dates.DateTime, firmness, TS::Vector{Dates.DateTime}, context::AbstractContext)
    @info("-----Entrée dans la fenêtre opérationnelle-----")
    return nothing #result
end
function compute_firmness(runnable::EnterFO,
                        ech::Dates.DateTime, next_ech::Union{Nothing,Dates.DateTime},
                        TS::Vector{Dates.DateTime}, context::AbstractContext)
    return Firmness()
end

