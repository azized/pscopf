using Dates
using Printf

using .Networks

###############################################
# Sequence Generation
###############################################

struct SequenceGenerator <: AbstractDataGenerator
    network::Networks.Network #Not used for now (but potentially we can have specific operations at DMO horizons)
    target_timepoints::Vector{Dates.DateTime}
    horizon_timepoints::Vector{Dates.DateTime}
    management_mode::ManagementMode
end

function launch(seq_generator::SequenceGenerator)
    if seq_generator.management_mode == PSCOPF_MODE_1
        return gen_seq_mode1(seq_generator)

    elseif seq_generator.management_mode == PSCOPF_MODE_2
        return gen_seq_mode2(seq_generator)

    elseif seq_generator.management_mode == PSCOPF_MODE_3
        return gen_seq_mode3(seq_generator)
    end

    error("unsuppported mode : ", seq_generator.management_mode)
end

"""
    generate_sequence
"""
function generate_sequence(network::Networks.Network, target_timepoints::Vector{Dates.DateTime},
                            horizon_timepoints::Vector{Dates.DateTime}, management_mode::ManagementMode)
    generator = SequenceGenerator(network, target_timepoints, horizon_timepoints, management_mode)
    return launch(generator)
end


function gen_seq_mode1(seq_generator::SequenceGenerator)
    sequence = Sequence()
    first_ts = seq_generator.target_timepoints[1]
    fo_startpoint = first_ts - get_fo_length(seq_generator.management_mode)

    for ech in seq_generator.horizon_timepoints
        if ech <  first_ts
            if ech < fo_startpoint
                add_step!(sequence, EnergyMarket, ech)
                add_step!(sequence, TSOOutFO, ech)

            elseif ech == fo_startpoint
                add_step!(sequence, EnergyMarketAtFO, ech)
                add_step!(sequence, EnterFO, ech)
                add_step!(sequence, TSOInFO, ech)

            else
                add_step!(sequence, TSOInFO, ech)
            end
        elseif first_ts < ech
            msg = @sprintf(("Error when generating sequence: ech (%s) is after target timepoint (%s)."),
                            ech, first_ts)
            throw( error(msg) )
        end
    end

    add_step!(sequence, Assessment, seq_generator.horizon_timepoints[end])

    return sequence
end


function gen_seq_mode2(seq_generator::SequenceGenerator)
    sequence = Sequence()
    first_ts = seq_generator.target_timepoints[1]
    fo_startpoint = first_ts - get_fo_length(seq_generator.management_mode)

    preceding_ech = nothing
    for ech in seq_generator.horizon_timepoints
        if ech <  first_ts
            if ech < fo_startpoint
                add_step!(sequence, EnergyMarket, ech)
                add_step!(sequence, TSOOutFO, ech)

            elseif ech == fo_startpoint
                add_step!(sequence, EnergyMarketAtFO, ech)
                add_step!(sequence, EnterFO, ech)
                add_step!(sequence, TSOOutFO, ech)

            elseif preceding_ech == fo_startpoint
                add_step!(sequence, EnergyMarket(EnergyMarketConfigs(REF_SCHEDULE_TYPE=TSO())), ech)
                add_step!(sequence, TSOBilevel(TSOBilevelConfigs(REF_SCHEDULE_TYPE_IN_TSO=TSO())), ech) #Ref can be market cause it respects the preceding TSO
                add_step!(sequence, BalanceMarket(EnergyMarketConfigs(REF_SCHEDULE_TYPE=TSO())), ech)
            else
                add_step!(sequence, EnergyMarket(EnergyMarketConfigs(REF_SCHEDULE_TYPE=Market())), ech)
                # NOTE: The TSOBilevel considers the EnergyMarket both as a reference for deltas and for decided schedules
                #we can keep this cause the EnergyMarket respects the decisions of the previous BalanceMarket
                add_step!(sequence, TSOBilevel(TSOBilevelConfigs(REF_SCHEDULE_TYPE_IN_TSO=Market())), ech)
                add_step!(sequence, BalanceMarket(EnergyMarketConfigs(REF_SCHEDULE_TYPE=Market())), ech)
            end
        elseif first_ts < ech
            msg = @sprintf(("Error when generating sequence: ech (%s) is after target timepoint (%s)."),
                            ech, first_ts)
            throw( error(msg) )
        end
        preceding_ech = ech
    end

    add_step!(sequence, Assessment, seq_generator.horizon_timepoints[end])

    return sequence
end


function gen_seq_mode3(seq_generator::SequenceGenerator)
    sequence = Sequence()
    first_ts = seq_generator.target_timepoints[1]
    fo_startpoint = first_ts - get_fo_length(seq_generator.management_mode)

    for ech in seq_generator.horizon_timepoints
        if ech <  first_ts
            if ech < fo_startpoint
                add_step!(sequence, EnergyMarket, ech)
                add_step!(sequence, TSOOutFO, ech)

            elseif ech == fo_startpoint
                add_step!(sequence, EnergyMarket, ech)
                add_step!(sequence, EnterFO, ech)
                add_step!(sequence, TSOBilevel(TSOBilevelConfigs(REF_SCHEDULE_TYPE_IN_TSO=TSO())), ech)

            else
                #TODO : check if this is an EnergyMarket or BalanceMarket or another implem
                add_step!(sequence, BalanceMarket, ech)
            end
        elseif first_ts < ech
            msg = @sprintf(("Error when generating sequence: ech (%s) is after target timepoint (%s)."),
                            ech, first_ts)
            throw( error(msg) )
        end
    end

    add_step!(sequence, Assessment, seq_generator.horizon_timepoints[end])

    return sequence
end
