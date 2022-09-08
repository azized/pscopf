
module PSCOPFio

using ..PSCOPF
using ..Networks

using Dates
using Printf
using DataStructures

OUTPUT_PREFIX = "pscopf_out_"

##########################
#   Readers
##########################


function read_buses_from_ptdf(data_folder::String)::Set{String}
    buses_ids = Set{String}();
    open(joinpath(data_folder, "pscopf_ptdf.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);
                push!(buses_ids, buffer[3])
            end
        end
    end
    return buses_ids
end

function read_buses_from_units(data_folder::String)::Set{String}
    buses_ids = Set{String}();
    open(joinpath(data_folder, "pscopf_units.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);
                push!(buses_ids, buffer[3])
            end
        end
    end
    return buses_ids
end

function read_buses!(network::Network, data::String)
    buses_ids::Set{String} = union(
        read_buses_from_ptdf(data),
        read_buses_from_units(data)
    )
    Networks.add_new_buses!(network, buses_ids);
end


function read_branches_from_ptdf(data_folder::String, default_limit=0.)::Dict{String, Networks.Limits}
    branches = Dict{String, Networks.Limits}();
    open(joinpath(data_folder, "pscopf_ptdf.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);
                network_case = buffer[1]
                branch_id = buffer[2]
                limits_l = get!(branches, branch_id, Networks.Limits())
                push!(limits_l, network_case => default_limit)
            end
        end
    end
    return branches
end

function read_branches_from_limits(data_folder::String)::Dict{String, Networks.Limits}
    branches = Dict{String, Networks.Limits}();
    open(joinpath(data_folder, "pscopf_limits.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);
                network_case = buffer[1]
                branch_id = buffer[2]
                limit = parse(Float64, buffer[3]);
                limits_l = get!(branches, branch_id, Networks.Limits())
                push!(limits_l, network_case => limit)
            end
        end
    end
    return branches
end

function read_branches!(network::Network, data::String)
    default_limit = 0.
    # merge order is important to prioritise limits data!
    branches::Dict{String, Networks.Limits} = merge(
        read_branches_from_ptdf(data, default_limit),
        read_branches_from_limits(data)
    )

    for (id,limits) in branches
        Networks.add_new_branch!(network, id, limits);
    end
end

function read_ptdf!(network::Network, data::String, filename="pscopf_ptdf.txt")
    open(joinpath(data, filename), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);
                ptdf_case = buffer[1]
                if (PSCOPF.get_config("CONSIDER_N_1") || ptdf_case==Networks.BASECASE)
                    branch_id = buffer[2]
                    bus_id = buffer[3]
                    ptdf_value = parse(Float64, buffer[4])
                    Networks.add_ptdf_elt!(network, branch_id, bus_id, ptdf_value, ptdf_case)
                end
            end
        end
    end
end

function read_generators!(network, data)
    open(joinpath(data, "pscopf_units.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);

                generator_id = buffer[1]
                gen_type = parse(Networks.GeneratorType, buffer[2])
                gen_bus_id = buffer[3]
                pmin = parse(Float64, buffer[4])
                pmax = parse(Float64, buffer[5])
                start_cost = parse(Float64, buffer[6])
                prop_cost = parse(Float64, buffer[7])
                dmo = Dates.Second(parse(Float64, buffer[8]))
                dp = Dates.Second(parse(Float64, buffer[9]))

                Networks.add_new_generator_to_bus!(network, gen_bus_id,
                                        generator_id, gen_type, pmin, pmax, start_cost, prop_cost, dmo, dp)
            end
        end
    end
end

function read_uncertainties_distributions(data)
    result = PSCOPF.UncertaintiesDistribution()
    open(joinpath(data, "uncertainties_distribution.txt"), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln);

                id = buffer[1]
                min_value = parse(Float64, buffer[2])
                max_value = parse(Float64, buffer[3])
                mu = parse(Float64, buffer[4])
                sigma = parse(Float64, buffer[5])
                time_factor = parse(Float64, buffer[6])

                PSCOPF.add_uncertainty_distribution!(result,
                                        id, min_value, max_value, mu, sigma,
                                        time_factor)
            end
        end
    end

    return result
end

function read_uncertainties(data, filename="pscopf_uncertainties.txt")
    result = PSCOPF.Uncertainties()
    open(joinpath(data, filename), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln)
                # "name", "ts", "ech", "scenario", "value"))
                name = buffer[1]
                ts = Dates.DateTime(buffer[2])
                ech = Dates.DateTime(buffer[3])
                scenario = buffer[4]
                value = parse(Float64, buffer[5])
                PSCOPF.add_uncertainty!(result, ech, name, ts, scenario, value)
            end
        end
    end
    return result
end

function read_initial_state(data, filename="pscopf_init.txt")
    result = SortedDict{String, PSCOPF.GeneratorState}()
    open(joinpath(data, filename), "r") do file
        for ln in eachline(file)
            # don't read commentted line
            if ln[1] != '#'
                buffer = PSCOPF.split_with_space(ln)
                gen_id = buffer[1]
                state = parse(PSCOPF.GeneratorState, buffer[2])
                result[gen_id] = state
            end
        end
    end
    return result
end

##########################
#   Writers
##########################

function write(dir_path::String, network::Networks.Network;
            ignore_units=false, ignore_limits=false, ignore_ptdf=false)
    mkpath(dir_path)
    # if !isdir(dir_path)
    #     mkpath(dir_path)
    # else
    #     msg = @sprintf("data folder `%s` already exists!", dir_path)
    #     error(msg)
    # end
    #units
    if !ignore_units
        write(dir_path, network.generators)
    end
    #limits
    if !ignore_limits
        write(dir_path, network.branches)
    end
    #ptdf
    if !ignore_ptdf
        write(dir_path, network.ptdf)
    end
end

function write(dir_path::String, generators::SortedDict{String, Networks.Generator})
    output_file_l = joinpath(dir_path, "pscopf_units.txt")
    open(output_file_l, "w") do file_l
        Base.write(file_l, @sprintf("#%24s%16s%25s%16s%16s%16s%16s%16s%16s\n",
                    "name", "type", "bus_id", "minP","maxP", "start", "prop", "dmo(s)", "dp(s)"))
        for (id_l, generator_l) in generators
            Base.write(file_l, @sprintf("%25s%16s%25s%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n",
                                    Networks.get_id(generator_l),
                                    Networks.get_type(generator_l),
                                    Networks.get_bus_id(generator_l),
                                    Networks.get_p_min(generator_l),
                                    Networks.get_p_max(generator_l),
                                    Networks.get_start_cost(generator_l),
                                    Networks.get_prop_cost(generator_l),
                                    Dates.value(Networks.get_dmo(generator_l)),
                                    Dates.value(Networks.get_dp(generator_l))
                                    )
                    )
        end
    end
end

function write(dir_path::String, branches::SortedDict{String, Networks.Branch})
    output_file_l = joinpath(dir_path, "pscopf_limits.txt")
    open(output_file_l, "w") do file_l
        Base.write(file_l, @sprintf("#%24s %25s %16s\n", "network_case", "branch", "limit"))
        for (_, branch_l) in branches
            for (network_case_l, limit_l) in Networks.get_limit(branch_l)
                Base.write(file_l, @sprintf("%25s %25s %16.8E\n",
                                        network_case_l,
                                        Networks.get_id(branch_l),
                                        limit_l
                                        )
                        )
            end
        end
    end
end


function write(dir_path::String, ptdf_dict::Networks.PTDFDict)
    output_file_l = joinpath(dir_path, "pscopf_ptdf.txt")
    open(output_file_l, "w") do file_l
        Base.write(file_l, @sprintf("#%24s %16s\n", "REF_BUS", "unknown"))
        Base.write(file_l, @sprintf("#%24s %25s %25s %16s\n", "case", "branch", "bus", "value"))
        for (case_l, ptdf) in ptdf_dict
            for (branch_id_l, _) in ptdf
                for (bus_id_l, val_l) in ptdf[branch_id_l]
                    Base.write(file_l, @sprintf("%24s %25s %25s %16.8E\n",
                                            case_l,
                                            branch_id_l,
                                            bus_id_l,
                                            val_l
                                            )
                                )
                end
            end
        end
    end
end

function write(dir_path::String, gen_init::SortedDict{String, PSCOPF.GeneratorState})
    mkpath(dir_path)
    output_file_l = joinpath(dir_path, "pscopf_init.txt")
    open(output_file_l, "w") do file_l
        Base.write(file_l, @sprintf("#%24s%10s\n", "name", "state"))
        for (gen_id, state) in gen_init
            Base.write(file_l, @sprintf("%25s%10s\n", gen_id, state))
        end
    end
end

function write(dir_path::String, uncertainties::PSCOPF.Uncertainties)
    mkpath(dir_path)
    output_file_l = joinpath(dir_path, "pscopf_uncertainties.txt")
    open(output_file_l, "w") do file_l

        Base.write(file_l, @sprintf("#%24s%20s%20s%10s%16s\n", "name", "ts", "ech", "scenario", "value"))
        for (ech, _) in uncertainties
            for (nodal_injection_name, _) in uncertainties[ech]
                for (ts, _) in uncertainties[ech][nodal_injection_name]
                    for (scenario, value_l) in uncertainties[ech][nodal_injection_name][ts]
                        Base.write(file_l, @sprintf("%25s%20s%20s%10s%16.8E\n",
                                        nodal_injection_name,
                                        ts,
                                        ech,
                                        scenario,
                                        value_l
                                        )
                            )
                    end
                end
            end
        end

    end
end


function write_commitment_schedule(dir_path::String, schedule::PSCOPF.Schedule, prefix="")
    ech = schedule.decision_time
    commitment_filename_l = joinpath(dir_path, OUTPUT_PREFIX*prefix*"commitment_schedule.txt")
    open(commitment_filename_l, "a") do commitment_file_l
        if filesize(commitment_file_l) == 0
            Base.write(commitment_file_l, @sprintf("#%19s%10s%25s%20s%10s%6s%10s\n", "ech", "decider", "name", "ts", "scenario", "value", "firmness"))
        end
        for (gen_id, gen_schedule) in schedule.generator_schedules
            for (ts, uncertain_value) in gen_schedule.commitment
                firmness = PSCOPF.is_definitive(uncertain_value) ? "FIRM" : "FREE"
                for (scenario, value_l) in uncertain_value.anticipated_value
                    value_l = ismissing(value_l) ? -1. : value_l
                    Base.write(commitment_file_l, @sprintf("%20s%10s%25s%20s%10s%6s%10s\n",
                                        ech,
                                        schedule.decider_type,
                                        gen_id,
                                        ts,
                                        scenario,
                                        value_l,
                                        firmness
                                        )
                            )
                end
            end
        end
    end
end

function write_production_schedule(dir_path::String, schedule::PSCOPF.Schedule, prefix="")
    ech = schedule.decision_time
    schedule_filename_l = joinpath(dir_path, OUTPUT_PREFIX*prefix*"schedule.txt")
    open(schedule_filename_l, "a") do schedule_file_l
        if filesize(schedule_file_l) == 0
            Base.write(schedule_file_l, @sprintf("#%19s%10s%25s%20s%10s%16s%10s\n", "ech", "decider", "name", "ts", "scenario", "value", "firmness"))
        end
        for (gen_id, gen_schedule) in schedule.generator_schedules
            for (ts, uncertain_value) in gen_schedule.production
                firmness = PSCOPF.is_definitive(uncertain_value) ? "FIRM" : "FREE"
                for (scenario, value_l) in uncertain_value.anticipated_value
                    value_l = ismissing(value_l) ? -1. : value_l
                    Base.write(schedule_file_l, @sprintf("%20s%10s%25s%20s%10s%16.8E%10s\n",
                                        ech,
                                        schedule.decider_type,
                                        gen_id,
                                        ts,
                                        scenario,
                                        value_l,
                                        firmness
                                        )
                            )
                end
            end
        end
    end
end


function write_flows(dir_path::String, context::PSCOPF.AbstractContext, schedule::PSCOPF.Schedule, prefix="")
    ech = schedule.decision_time
    flows_filename_l = joinpath(dir_path, OUTPUT_PREFIX*prefix*"flows.txt")

    flows = PSCOPF.compute_flows(context, schedule)
    open(flows_filename_l, "a") do flows_file_l
        if filesize(flows_file_l) == 0
            Base.write(flows_file_l, @sprintf("#%19s%10s%25s%20s%10s%16s\n", "ech", "decider", "branch_name", "ts", "scenario", "value"))
        end
        for ((branch_id, ts, scenario), flow_value) in flows
            Base.write(flows_file_l, @sprintf("%20s%10s%25s%20s%10s%16.8E\n",
                                ech,
                                schedule.decider_type,
                                branch_id,
                                ts,
                                scenario,
                                flow_value
                                )
                    )
        end
    end
end

function write(context::PSCOPF.AbstractContext, schedule::PSCOPF.Schedule, prefix="")
    dir_path = context.out_dir
    if !isnothing(dir_path)
        mkpath(dir_path)
        write_commitment_schedule(dir_path, schedule, prefix)
        write_production_schedule(dir_path, schedule, prefix)
        write_flows(dir_path, context, schedule, prefix)
    end
end

function write_full(context::PSCOPF.AbstractContext, schedule::PSCOPF.Schedule, prefix="")
    dir_path = context.out_dir
    if !isnothing(dir_path)
        mkpath(dir_path)
        _write(dir_path, schedule, context, prefix)
    end
end

function _write(dir_path::String, schedule::PSCOPF.Schedule, context::PSCOPF.AbstractContext, prefix="")
    ech = schedule.decision_time
    decider = schedule.decider_type
    schedule_filename_l = joinpath(dir_path, OUTPUT_PREFIX*prefix*"full_schedule.txt")
    open(schedule_filename_l, "a") do schedule_file_l
        if filesize(schedule_file_l) == 0
            Base.write(schedule_file_l, @sprintf("#%19s%10s%25s%10s%20s%10s%16s%16s%16s%16s%15s%8s%15s\n",
                                            "ech", "decider", "unit/bus", "PROD/LOAD", "ts", "scenario",
                                            "inject_value", "power_value", "load_shed", "power_capped", "DP_firmness",
                                            "ON/OFF", "DMO_firmness",
                                            ))
        end

        TS =  PSCOPF.get_target_timepoints(context)
        scenarios = PSCOPF.get_scenarios(context)
        uncertainties_at_ech = PSCOPF.get_uncertainties(context, ech)
        network = PSCOPF.get_network(context)

        for generator in Networks.get_generators(network)
            gen_id = Networks.get_id(generator)
            gen_schedule = PSCOPF.get_sub_schedule(schedule, gen_id)

            load_shed = 0. #this is for buses only
            for ts in TS
                power_uncertain_value = PSCOPF.get_prod_uncertain_value(gen_schedule, ts)
                power_firmness = PSCOPF.is_definitive(power_uncertain_value) ? "FIRM" : "FREE"
                if Networks.needs_commitment(generator)
                    commitment_uncertain_value = PSCOPF.get_commitment_uncertain_value(gen_schedule, ts)
                    commitment_firmness = PSCOPF.is_definitive(commitment_uncertain_value) ? "FIRM" : "FREE"
                else
                    commitment_firmness = "X"
                end
                for scenario in scenarios
                    prod_value = PSCOPF.get_prod_value(gen_schedule, ts, scenario)
                    if Networks.is_limitable(generator)
                        power_value = PSCOPF.get_uncertainties(uncertainties_at_ech, gen_id, ts, scenario)
                        power_capped = PSCOPF.get_capping(schedule, gen_id, ts, scenario)
                        inject_value = prod_value - power_capped #TODO : add this to schedule and compute it by each Model
                                                                 #FIXME : for TSOModel the prod_value take account of capping already
                    else
                        inject_value = power_value = prod_value
                        power_capped = 0.
                    end
                    if Networks.needs_commitment(generator)
                        commitment_value = PSCOPF.get_commitment_value(gen_schedule, ts, scenario)
                    else
                        commitment_value = "X"
                    end

                    Base.write(schedule_file_l, @sprintf("%20s%8s%25s%10s%20s%10s%16.8E%16.8E%16.8E%16.8E%15s%8s%15s\n",
                                                    ech, decider, gen_id, "PROD", ts, scenario,
                                                    inject_value, power_value, load_shed, power_capped, power_firmness,
                                                    commitment_value, commitment_firmness)
                                )
                end
            end
        end

        for bus_id in Networks.get_id.(Networks.get_buses(network))
            for ts in TS
                for scenario in scenarios
                    power_value = PSCOPF.get_uncertainties(uncertainties_at_ech, bus_id, ts, scenario)
                    load_shed = PSCOPF.get_loss_of_load(schedule, bus_id, ts, scenario)
                    inject_value = power_value - load_shed

                    Base.write(schedule_file_l, @sprintf("%20s%8s%25s%8s%20s%10s%16.8E%16.8E%16.8E%16.8E%8s%6s%8s\n",
                                                        ech, decider, bus_id, "LOAD", ts, scenario,
                                                        inject_value, power_value, load_shed, 0., "X",
                                                        "X", "X")
                        )
                end
            end
        end
    end
end

end #module PSCOPFio
