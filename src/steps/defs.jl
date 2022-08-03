

struct OptimResult end

struct TSO <: DeciderType end
struct Market <: DeciderType end
struct Utilitary <: DeciderType end
struct Assess <: DeciderType end
DeciderType(::Type{<:AbstractRunnable}) = Utilitary()
DeciderType(::Type{<:AbstractTSO}) = TSO()
DeciderType(::Type{<:AbstractMarket}) = Market()
DeciderType(::AbstractRunnable) = Utilitary()
DeciderType(::AbstractAssessment) = Assess()
DeciderType(::AbstractTSO) = TSO()
DeciderType(::AbstractMarket) = Market()

is_tso(x::T) where {T} = is_tso(DeciderType(T))
is_tso(::DeciderType) = false
is_tso(::TSO) = true

is_market(x::T) where {T} = is_market(DeciderType(T))
is_market(::DeciderType) = false
is_market(::Market) = true

is_utilitary(x::T) where {T} = is_utilitary(DeciderType(T))
is_utilitary(::DeciderType) = false
is_utilitary(::Utilitary) = true

function Base.string(decider_type::DeciderType)
    if is_market(decider_type)
        return "market"
    elseif is_tso(decider_type)
        return "tso"
    elseif typeof(decider_type) <: Utilitary
        return "util"
    else
        throw( error("Undefined conversion of DeciderType `", decider_type, "` to a string") )
    end
end
