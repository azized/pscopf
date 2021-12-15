

export get_ts_s;
export get_units_by_kind;
export get_units_by_bus;
export get_bus;
export sc_opf;
export print_nz;

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
                throw(e_xpress)
            end
        end
    else
        throw(e_xpress)
    end
end
println("optimizer: ", OPTIMIZER)

using Dates
import Statistics

function get_bus(launcher::Launcher, names::Set{String})
    result = Set{String}();
    for name in names
        if !haskey(launcher.gen_type_bus, name)
            push!(result, name);
        end
    end
    return result;
end

function get_sorted_ech(launcher::Launcher)
    ech_set = Set{DateTime}();
    for uncertainty in launcher.uncertainties
        ech = uncertainty[1][U_ECH];
        push!(ech_set, ech);
    end
    return sort(collect(ech_set));
end

function get_ts_s_name(launcher::Launcher, ech_p::DateTime)
    ts_set = Set{DateTime}();
    s_set =  Set{String}();
    name_set = Set{String}();
    for uncertainty in launcher.uncertainties
        ech_l = uncertainty[1][U_ECH];
        if ech_l == ech_p
            ts = uncertainty[1][U_H];
            s = uncertainty[1][U_SCENARIO];
            name = uncertainty[1][U_NAME];
            push!(ts_set, ts);
            push!(s_set, s);
            push!(name_set, name);
        end
    end
    return sort(collect(ts_set)), sort(collect(s_set)), name_set;
end

function get_ech_ts_s_name(launcher::Launcher)
    ech_set = Set{DateTime}();
    ts_set = Set{DateTime}();
    s_set =  Set{String}();
    name_set = Set{String}();
    for uncertainty in launcher.uncertainties
        ts = uncertainty[1][U_H];
        ech = uncertainty[1][U_ECH];
        s = uncertainty[1][U_SCENARIO];
        name = uncertainty[1][U_NAME];
        push!(ech_set, ech);
        push!(ts_set, ts);
        push!(s_set, s);
        push!(name_set, name);
    end
    return sort(collect(ech_set)), sort(collect(ts_set)), sort(collect(s_set)), name_set;
end

function get_units_by_kind(launcher::Launcher)
    result = Dict{String,Dict{String,String}}();
    result[K_LIMITABLE] = Dict{String,String}();
    result[K_IMPOSABLE] = Dict{String,String}();
    for gen in launcher.gen_type_bus
        kind = gen[2][1];
        push!(result[kind],  gen[1] => gen[2][2]);
    end 
    return result;
end

function get_units_by_bus(launcher::Launcher, buses::Set{String})
    result = Dict{String,Dict{String,Vector{String}}}();
    result[K_LIMITABLE] = Dict{String,Vector{String}}([bus => Vector{String}() for bus in buses]);
    result[K_IMPOSABLE] = Dict{String,Vector{String}}([bus => Vector{String}() for bus in buses]);
    for gen in launcher.gen_type_bus
        name = gen[1];
        kind = gen[2][1];
        bus  = gen[2][2];
        tmp = get(result[kind], bus, Vector{String}())
        push!(tmp, name);
        result[kind][bus] = tmp;
    end 
    return result;
end
function add_limitable!(launcher::Launcher, ech::DateTime, model, units_by_kind, TS, S)
    p_lim = Dict{Tuple{String,DateTime},VariableRef}();
    is_limited = Dict{Tuple{String,DateTime,String},VariableRef}();
    p_enr = Dict{Tuple{String,DateTime,String},VariableRef}();
    is_limited_x_p_lim = Dict{Tuple{String,DateTime,String},VariableRef}();
    c_lim = Dict{Tuple{String,DateTime,String},VariableRef}();
    if ! launcher.NO_LIMITABLE
        for kvp in units_by_kind[K_LIMITABLE], ts in TS
            gen = kvp[1];
            pLimMax = 0;
            p_lim[gen, ts] =  @variable(model, base_name = @sprintf("p_lim[%s,%s]", gen, ts), lower_bound = 0);
            for s in S
                pLimMax = max(pLimMax, launcher.uncertainties[gen, s, ts, ech]);
            end
            for s in S
                name =  @sprintf("is_limited[%s,%s,%s]", gen, ts, s);
                is_limited[gen, ts, s] = @variable(model, base_name = name, binary = true);

                name =  @sprintf("is_limited_x_p_lim[%s,%s,%s]", gen, ts, s);
                is_limited_x_p_lim[gen, ts, s] = @variable(model, base_name = name, lower_bound = 0);

                name =  @sprintf("c_lim[%s,%s,%s]", gen, ts, s);
                c_lim[gen, ts, s] = @variable(model, base_name = name, lower_bound = 0);

                @constraint(model, is_limited_x_p_lim[gen, ts, s] <= p_lim[gen, ts]);
                @constraint(model, is_limited_x_p_lim[gen, ts, s] <= is_limited[gen, ts, s] * pLimMax);
                @constraint(model, is_limited_x_p_lim[gen, ts, s] + pLimMax * (1 - is_limited[gen, ts, s]) - p_lim[gen, ts] >= 0);

                p0 = launcher.uncertainties[gen, s, ts, ech];
                name =  @sprintf("p_enr[%s,%s,%s]", gen, ts, s);
                p_enr[gen, ts, s] = @variable(model, base_name = name, lower_bound = 0, upper_bound = p0);
                @constraint(model, p_enr[gen, ts, s] == (1 - is_limited[gen, ts, s]) * p0 + is_limited_x_p_lim[gen, ts, s]);
                @constraint(model, p_enr[gen, ts, s] <= p_lim[gen, ts]);

                @constraint(model, c_lim[gen, ts, s] >= p0 - p_lim[gen, ts]);

                if launcher.NO_LIMITATION
                    @constraint(model, is_limited[gen, ts, s] == 0);
                end
            end
        end
    end
    return Workflow.LimitableModeler(p_lim, is_limited, p_enr, is_limited_x_p_lim, c_lim);
    # return p_lim, is_limited, is_limited_x_p_lim, c_lim;
end

"""
    is_already_fixed(ech_p, ts_p, dmo_p)

returns true if the production level for time step ts_p for a unit having a delay equal to dmo_p can no longer be changed by the time ech_p

# Arguments
- `ech_p::DateTime` : is the time we decide (the optimisation launch time)
- `ts_p::DateTime` : is the production time
- `dmo_p` : is the necessary time (in seconds) for the unit to start producing. Should be losslessly convertible to Int64.
"""
function is_already_fixed(ech_p::DateTime, ts_p::DateTime, dmo_p)
    return (ts_p - Dates.Second(dmo_p)) < ech_p
end

"""
    is_to_decide(ech_p, ts_p, dmo_p)

returns true if the production level for time step ts_p for a unit having a delay equal to dmo_p must definately be decided at time ech_p

# Arguments
- `ech_p::DateTime` : is the time we decide (the optimisation launch time)
- `ts_p::DateTime` : is the production time
- `dmo_p` : is the necessary time (in seconds) for the unit to start producing. Should be losslessly convertible to Int64.
"""
function is_to_decide(ech_p::DateTime, ts_p::DateTime, dmo_p)
    return (ts_p - Dates.Second(dmo_p)) == ech_p
end

function add_imposable!(launcher::Launcher, ech, model,  units_by_kind, TS, S)
    p_imposable = Dict{Tuple{String,DateTime,String},VariableRef}();
    p_is_imp = Dict{Tuple{String,DateTime,String},VariableRef}();
    p_is_imp_and_on = Dict{Tuple{String,DateTime,String},VariableRef}();
    p_imp = Dict{Tuple{String,DateTime,String},VariableRef}();
    p_start = Dict{Tuple{String,DateTime,String},VariableRef}();
    p_on = Dict{Tuple{String,DateTime,String},VariableRef}();
    c_imp_pos  = Dict{Tuple{String,DateTime,String},VariableRef}();
    c_imp_neg   = Dict{Tuple{String,DateTime,String},VariableRef}();

    if ! launcher.NO_IMPOSABLE
        for kvp in units_by_kind[K_IMPOSABLE]
            gen = kvp[1];
            pMin = launcher.units[gen][1];
            pMax = launcher.units[gen][2];
            dmo_l = launcher.units[gen][5];
            for ts in TS
                for s in S
                    name =  @sprintf("p_imp[%s,%s,%s]", gen, ts, s);
                    p_imp[gen, ts, s] = @variable(model, base_name = name);

                    name =  @sprintf("p_is_imp[%s,%s,%s]", gen, ts, s);
                    p_is_imp[gen, ts, s] = @variable(model, base_name = name, binary = true);

                    name =  @sprintf("p_start[%s,%s,%s]", gen, ts, s);
                    p_start[gen, ts, s] = @variable(model, base_name = name, binary = true);

                    name =  @sprintf("p_on[%s,%s,%s]", gen, ts, s);
                    p_on[gen, ts, s] = @variable(model, base_name = name, binary = true);

                    name =  @sprintf("p_is_imp_and_on[%s,%s,%s]", gen, ts, s);
                    p_is_imp_and_on[gen, ts, s] = @variable(model, base_name = name, binary = true);

                    # if in(gen, ["park_city", "alta", "sundance"])
                    #     println("WARNING");
                    #     @constraint(model,  p_on[gen, ts]==0);
                    #     @constraint(model,  p_is_imp[gen, ts]==0);
                    # end

                    p0 = launcher.uncertainties[gen, s, ts, ech];
                    name =  @sprintf("p_imposable[%s,%s,%s]", gen, ts, s);
                    p_imposable[gen, ts, s] = @variable(model, base_name = name, lower_bound = 0);

                    name =  @sprintf("c_imp_pos[%s,%s,%s]", gen, ts, s);
                    c_imp_pos[gen, ts, s] = @variable(model, base_name = name, lower_bound = 0);

                    name =  @sprintf("c_imp_neg[%s,%s,%s]", gen, ts, s);
                    c_imp_neg[gen, ts, s] = @variable(model, base_name = name, lower_bound = 0);

                    @constraint(model, p_imposable[gen, ts, s] == p0 + c_imp_pos[gen, ts, s] - c_imp_neg[gen, ts, s]);

                    @constraint(model, p_imposable[gen, ts, s] == (1 - p_is_imp[gen, ts, s]) * p0 + p_imp[gen, ts, s]);
                    @constraint(model, p_imp[gen, ts, s] <= pMax * p_is_imp_and_on[gen, ts, s]);
                    @constraint(model, p_imp[gen, ts, s] >= pMin * p_is_imp_and_on[gen, ts, s]);

                    @constraint(model, p_is_imp_and_on[gen, ts, s] <= p_is_imp[gen, ts, s]);
                    @constraint(model, p_is_imp_and_on[gen, ts, s] <= p_on[gen, ts, s]);
                    @constraint(model, 1 + p_is_imp_and_on[gen, ts, s] >= p_on[gen, ts, s] + p_is_imp[gen, ts, s]);

                    if is_already_fixed(ech, ts, dmo_l)
                    # it is too late to change the production level
                        #FIXME : Discuss uncertainties : cause if the uncertainties are different for this scenario it will be considered as imposed
                        @printf("%s: unit %s is already decided for timestep %s.\n", ech, gen, ts)
                        prev0 = launcher.previsions[gen, ts, ech];
                        @constraint(model, p_imposable[gen, ts, s] == prev0);
                    elseif is_to_decide(ech, ts, dmo_l)
                    # must decide now on production level
                        @printf("%s: unit %s must be fixed for timestep %s.\n", ech, gen, ts)
                        for ((_,_,s_other_l), _) in filter(x ->  ( (x[1][1]==gen) && (x[1][2]==ts) ), p_imposable)
                            #iterate on other scenarios for which we already defined a variable : for s' in S st s' < s
                            @constraint(model, p_imposable[gen, ts, s] == p_imposable[gen, ts, s_other_l]);
                        end
                    else
                        @printf("%s: still early to fix unit %s for timestep %s.\n", ech, gen, ts)
                    end
                end
            end
            for s in S
                for (i,ts) in enumerate(TS)
                    if i > 1
                        ts_1 = TS[i - 1];
                        @constraint(model, p_start[gen, ts, s] <= p_on[gen, ts, s]);
                        @constraint(model, p_start[gen, ts, s] <= 1 - p_on[gen, ts_1, s]);
                        @constraint(model, p_start[gen, ts, s] >= p_on[gen, ts, s] - p_on[gen, ts_1, s]);
                    else
                        prev0 = launcher.previsions[gen, ts, ech];
                        println(@sprintf("%s %s is at %f\n", gen, ts, prev0));
                        if abs(prev0) < 1
                            @constraint(model, p_start[gen, ts, s] == p_on[gen, ts, s]);
                        end
                    end
                end
            end
        end
    end
    return Workflow.ImposableModeler(p_imposable, p_is_imp, p_imp, p_is_imp_and_on, p_start, p_on, c_imp_pos, c_imp_neg);
end

function add_eod_constraint!(launcher::Launcher, ech, model, units_by_bus, TS, S, BUSES, v_lim::Workflow.LimitableModeler, v_imp::Workflow.ImposableModeler, v_res::Workflow.ReserveModeler, netloads )
    eod_expr = Dict([(ts, s) => AffExpr(0) for ts in TS, s in S]);
    for ts in TS, s in S
        for bus in BUSES
            if ! launcher.NO_IMPOSABLE
                for gen in units_by_bus[K_IMPOSABLE][bus]
                    eod_expr[ts, s] += v_imp.p_imposable[gen, ts, s];
                end
            end
            if ! launcher.NO_LIMITABLE
                for gen in units_by_bus[K_LIMITABLE][bus]
                    p0 = launcher.uncertainties[gen, s, ts, ech];
                    # println(gen, " ", ts, " ", s, " ", p0);
                    #eod_expr[ts, s] +=  ((1 -  v_lim.is_limited[gen, ts, s]) * p0 + v_lim.is_limited_x_p_lim[gen, ts, s]);
                    eod_expr[ts, s] += v_lim.p_enr[gen,ts,s]
                end
            end
        end
    end
    for ts in TS, s in S
        eod_expr[ts, s] += v_res.p_res_pos[ts, s];
        eod_expr[ts, s] -= v_res.p_res_neg[ts, s];
    end
    
    # println(launcher.uncertainties)
    for ts in TS, s in S
        load = 0;
        for bus in BUSES
            load += netloads[bus, ts, s];
        end
        @constraint(model, eod_expr[ts, s] == load);
    end
end

function add_reserve!(launcher::Launcher, ech, model, TS, S, p_res_min, p_res_max)
    p_res_pos = Dict{Tuple{DateTime,String},VariableRef}();
    p_res_neg = Dict{Tuple{DateTime,String},VariableRef}();
    for ts in TS, s in S
        name =  @sprintf("p_res_pos[%s,%s]", ts, s);
        p_res_pos[ts, s] = @variable(model, base_name = name, lower_bound = 0, upper_bound = p_res_max);
        name =  @sprintf("p_res_neg[%s,%s]", ts, s);
        p_res_neg[ts, s] = @variable(model, base_name = name, lower_bound = 0, upper_bound = -p_res_min);
    end
    return Workflow.ReserveModeler(p_res_pos, p_res_neg);
end

function add_obj!(launcher::Launcher, model, TS, S, v_lim::Workflow.LimitableModeler, v_imp::Workflow.ImposableModeler, v_res::Workflow.ReserveModeler)
    eod_obj_l = AffExpr(0);

    penalties_obj_l = AffExpr(0);
    lim_cost_obj_l = AffExpr(0);
    imp_prop_cost_obj_l = AffExpr(0);
    imp_starting_cost_obj_l = AffExpr(0);

    w_scenario_l = 1 / length(S);

    for x in v_lim.c_lim
        # implicit looping on #ts, #s and gen
        gen = x[1][1];
        cProp = launcher.units[gen][4];
        lim_cost_obj_l += w_scenario_l * cProp * x[2];
    end
    if length(v_lim.is_limited) > 0
        penalties_obj_l += sum(1e-4 * x[2] for x in v_lim.is_limited);
    end

    for x in v_imp.p_start
        gen = x[1][1];
        cStart = launcher.units[gen][3];
        imp_starting_cost_obj_l += w_scenario_l * cStart * x[2];
    end
    for x in v_imp.c_imp_pos
        gen = x[1][1];
        cProp = launcher.units[gen][4];
        imp_prop_cost_obj_l += w_scenario_l * cProp * x[2];
    end
    for x in v_imp.c_imp_neg
        gen = x[1][1];
        cProp = launcher.units[gen][4];
        imp_prop_cost_obj_l += w_scenario_l * cProp * x[2];
    end

    for x in v_res.p_res_pos
        penalties_obj_l += (1e-3)*w_scenario_l*x[2];
    end
    for x in v_res.p_res_neg
        penalties_obj_l += (1e-3)*w_scenario_l* x[2];
    end

    eod_obj_l += penalties_obj_l
    eod_obj_l += lim_cost_obj_l
    eod_obj_l += imp_prop_cost_obj_l
    eod_obj_l += imp_starting_cost_obj_l
    obj_l = @objective(model, Min, eod_obj_l);

    return Workflow.ObjectiveModeler( penalties_obj_l, lim_cost_obj_l, imp_prop_cost_obj_l, imp_starting_cost_obj_l, obj_l)
end

function add_flow!(launcher::Workflow.Launcher, model, TS, S,  units_by_bus, v_lim::Workflow.LimitableModeler, v_imp::Workflow.ImposableModeler, netloads)
    v_flow = Dict{Tuple{String,DateTime,String},VariableRef}();

    ptdf_by_branch = Dict{String,Dict{String,Float64}}() ;
    for ptdf in launcher.ptdf
        branch, bus = ptdf[1]
        if ! haskey(ptdf_by_branch, branch)
            ptdf_by_branch[branch] = Dict{String,Float64}();
        end
        ptdf_by_branch[branch][bus] = ptdf[2];
    end
    println("ptdf_by_branch : ", ptdf_by_branch);
    for ts in TS, s in S
        for limit in launcher.limits
            branch = limit[1];
            v = limit[2]
            name =  @sprintf("p_flow[%s, %s,%s]", branch, ts, s);
            v_flow[branch, ts, s] = @variable(model, base_name = name, lower_bound = -v, upper_bound = v);
            f_ref = 0;
            xpr = AffExpr(0);
            for (bus, ptdf) in ptdf_by_branch[branch]
                f_ref -= ptdf * netloads[bus, ts, s];
                if ! launcher.NO_IMPOSABLE
                    for gen in units_by_bus[K_IMPOSABLE][bus]
                        xpr += v_imp.p_imposable[gen, ts, s] * ptdf;
                    end
                end
                if ! launcher.NO_LIMITABLE
                    for gen in units_by_bus[K_LIMITABLE][bus]
                        xpr += v_lim.p_enr[gen, ts, s] * ptdf;
                    end
                end
            end
            @constraint(model, v_flow[branch, ts, s] == f_ref + xpr);
        end
    end
    return v_flow;
end


"""
    sc_opf(launcher::Launcher, ech::DateTime, p_res_min, p_res_max)

Launch a single optimization iteration

# Arguments
- `launcher::Launcher` : the optimization launcher containing the necessary data
- `ech_p::DateTime` : the optimization launch time
- `p_res_min` : The minimum allowed reserve level
- `p_res_max` : The maximum allowed reserve level
"""
function sc_opf(launcher::Launcher, ech::DateTime, p_res_min, p_res_max)
    ##############################################################
    ### optimisation modelling sets
    ##############################################################
    TS, S, NAMES = Workflow.get_ts_s_name(launcher, ech);
    BUSES =  Workflow.get_bus(launcher, NAMES);
    units_by_kind = Workflow.get_units_by_kind(launcher);
    units_by_bus =  Workflow.get_units_by_bus(launcher, BUSES);

    println("NAMES : ", NAMES);
    println("BUSES : ", BUSES);
    println("units_by_bus : ", units_by_bus);

    println("ech: ", ech);
    println("Number of scenario ", length(S));
    println("Number of time step is ", length(TS));
    K_IMPOSABLE = Workflow.K_IMPOSABLE;
    K_LIMITABLE = Workflow.K_LIMITABLE;

    netloads = Dict([(bus, ts, s) => launcher.uncertainties[bus, s, ts, ech] for bus in BUSES, ts in TS, s in S]);
    eod_slack = Dict([(ts,s) => 0.0 for s in S, ts in TS]);
    factor = Dict([name => 1 for name in NAMES]);
    for bus in BUSES
        factor[bus] = -1;
    end
    for ((name, s, ts, ech), v) in launcher.uncertainties
        eod_slack[ts,  s] += factor[name] * v;
    end
    println("eod_slack is $eod_slack");
    ##############################################################
    # to be in a function ...
    ##############################################################

    model = Model();
    v_lim = add_limitable!(launcher, ech, model, units_by_kind, TS, S);
    p_lim = v_lim.p_lim;
    is_limited = v_lim.is_limited;
    is_limited_x_p_lim = v_lim.is_limited_x_p_lim;
    c_lim = v_lim.c_lim;

    v_imp = add_imposable!(launcher, ech, model, units_by_kind, TS, S);

    p_imposable = v_imp.p_imposable;
    p_is_imp = v_imp.p_is_imp;
    p_imp = v_imp.p_imp;
    p_start = v_imp.p_start;
    p_on = v_imp.p_on;

    v_res = add_reserve!(launcher, ech, model, TS, S, p_res_min, p_res_max);

    add_eod_constraint!(launcher, ech, model, units_by_bus, TS, S, BUSES, v_lim, v_imp, v_res, netloads);

    v_flow = add_flow!(launcher, model, TS, S, units_by_bus,  v_lim, v_imp, netloads);

    # println(units_by_bus)
    # println(eod_expr)
    obj_l = add_obj!(launcher,  model, TS, S, v_lim, v_imp, v_res);

    # println(model)
    set_optimizer(model, OPTIMIZER);
    write_to_file(model, @sprintf("%s/model_%s.lp", launcher.dirpath, ech))
    optimize!(model);

    println("end of optim.")
    if termination_status(model) == INFEASIBLE
        error("Model is infeasible")
    end
    println(termination_status(model))
    println(objective_value(model))

    # print_nz(p_imposable);
    # print_nz(p_is_imp);
    # print_nz(p_start);
    # print_nz(p_on);
    # print_nz(is_limited);

    println("NO_IMPOSABLE   : ", launcher.NO_IMPOSABLE);
    println("NO_LIMITABLE   : ", launcher.NO_LIMITABLE);
    println("NO_LIMITATION  : ", launcher.NO_LIMITATION);
    println("BUSES          : ", BUSES);
    println("NAMES          : ", NAMES);
    println("TS             : ", TS);
    println("S              : ", S);
    
    @printf("%30s[%s,%s] %4s %10s %10s %10s\n", "gen", "ts", "s", "b_enr", "p_enr", "p_lim_enr", "p0");
    if ! launcher.NO_LIMITABLE
        # println(launcher.uncertainties)
        for bus in BUSES, ts in TS, s in S, gen in units_by_bus[K_LIMITABLE][bus]
            b_lim = value(is_limited[gen, ts, s]);
            if b_lim > 0.5
                @printf("%10s[%s,%s] LIMITED\n", gen, ts, s);
            end
            p0 = launcher.uncertainties[gen, s, ts, ech];
            b_enr = value(is_limited[gen, ts, s]);
            p_enr = value((1 - is_limited[gen, ts, s]) * p0 + is_limited_x_p_lim[gen, ts, s]);
            p_lim_enr = value(p_lim[gen, ts]);
            # println(value(is_limited_x_p_lim[gen, ts, s]) - b_enr * p_lim_enr)
            if b_enr > 0.5
                @printf("%10s[%s,%s] %4d %10.3f %10.3f %10.3f\n", gen, ts, s, b_enr, p_enr, p_lim_enr, p0);
            end
            # println(gen, " : ", value((1 - is_limited[gen, ts, s]) * p0 + is_limited_x_p_lim[gen, ts, s]))
        end
    end
    if ! launcher.NO_IMPOSABLE
        # println(launcher.uncertainties)
        for bus in BUSES, ts in TS, s in S, gen in units_by_bus[K_IMPOSABLE][bus]
            b_start = value(p_start[gen, ts, s]);
            b_imp = value(p_is_imp[gen, ts, s]);
            if b_imp > 0.5
                @printf("%10s[%s, %s] IMPOSED\n", gen, ts, s);
            end
            if b_start > 0.5
                @printf("%10s[%s, %s] STARTED\n", gen, ts, s);
            end
            b_on = value(p_on[gen, ts, s]);
            if b_on > 0.5
                @printf("%10s[%s, %s] %10.3f\n", gen, ts, s, value(v_imp.p_imposable[gen, ts, s]));
            end
        end
    end
    println("v_res.p_res_pos");
    print_nz(v_res.p_res_pos);
    println("v_res.p_res_neg");
    print_nz(v_res.p_res_neg);
    println("v_res.v_flow");
    print_nz(v_flow);
    write_output_csv(launcher, ech, TS, S, v_lim, v_imp);
    write_extra_output_csv(launcher, ech, v_res, v_flow);

    return model, v_lim, v_imp, v_res, v_flow;
end

"""
update_schedule!(launcher_p::Launcher, next_ech_p::DateTime, ech_p::DateTime, v_lim_p::Workflow.LimitableModeler, v_imp_p::Workflow.ImposableModeler)

    Launch the optimization for multiple launch dates

    # Arguments
    - `launcher_p::Launcher` : the optimization launcher containing the necessary data. Attribute previsions will be updated.
    - `next_ech_p::DateTime` : The next date at which optimization will be relaunched
    - `ech_p::DateTime` : The last date we launch optimization at. (The values we are using for the update)
    - `v_lim_p::LimitableModeler` : container for the decision values for limitable units at time ech_p
    - `v_imp_p::ImposableModeler` : container for the decision values for imposable units at time ech_p
"""
function update_schedule!(launcher_p::Launcher, next_ech_p::DateTime, ech_p::DateTime, v_lim_p::Workflow.LimitableModeler, v_imp_p::Workflow.ImposableModeler)
    # keep previsionnal planning for dates != next_ech_p
    filter!( x -> x[1][3] !=  next_ech_p, launcher_p.previsions)

    #NOTE: This supposes TS(next_ech_p) included in TS(ech_p)
    TS, S, _ = get_ts_s_name(launcher_p, ech_p)
    units_by_kind_l = get_units_by_kind(launcher_p)

    for ts_l in TS
        for (gen_limitable_l,_) in units_by_kind_l[K_LIMITABLE]
            values_l = []
            for s_l in S
                push!(values_l, value(v_lim_p.p_enr[gen_limitable_l, ts_l, s_l]))
            end
            launcher_p.previsions[gen_limitable_l,ts_l,next_ech_p] = Statistics.mean(values_l)
        end
        for (gen_imposable_l,_) in units_by_kind_l[K_IMPOSABLE]
            values_l = []
            for s_l in S
                push!(values_l, value(v_imp_p.p_imposable[gen_imposable_l, ts_l, s_l]))
            end
            launcher_p.previsions[gen_imposable_l,ts_l,next_ech_p] = Statistics.mean(values_l)
        end
    end

    return launcher_p.previsions
end

function print_nz(variables)    
    for x in variables
        if abs(value(x[2])) > 1e-6
            println(x, " = ", value(x[2]));
        end
    end
end    

function clear_output_files(launcher)
    for output_file_l in CLEARED_OUTPUT
        output_path_l = joinpath(launcher.dirpath, output_file_l)
        if isfile(output_path_l)
            rm(output_path_l)
            println("removed ", output_path_l)
        end
    end
end

function write_output_csv(launcher::Workflow.Launcher, ech, TS, S, v_lim::Workflow.LimitableModeler, v_imp::Workflow.ImposableModeler)
    open(joinpath(launcher.dirpath, LIMITATION_CSV), "a") do file
        if filesize(file) == 0
            write(file, @sprintf("%s;%s;%s;%s;%s;%s;%s;%s;\n", "ech", "gen", "TS","S", "is_lim", "p_lim", "p0", "prev0"));
        end
        for x in v_lim.p_enr
            key = x[1];
            gen = key[1];
            ts =  key[2];
            s = key[3];

            is_lim = false;
            if value(v_lim.is_limited[key]) > 0.5
                is_lim = true;
            end
            p_sol = value(v_lim.p_lim[gen, ts]);
            p0 = launcher.uncertainties[gen, s, ts, ech];
            prev0 = launcher.previsions[gen, ts, ech];
            write(file, @sprintf("%s;%s;%s;%s;%s;%f;%f;%f\n", ech, gen, ts, s, is_lim, p_sol, p0, prev0));
        end
    end
    open(joinpath(launcher.dirpath, IMPOSITION_CSV), "a") do file
        if filesize(file) == 0
            write(file, @sprintf("%s;%s;%s;%s;%s;%s;%s;%s;\n", "ech", "gen", "TS","S", "is_imp", "p_imp", "p0", "prev0"));
        end
        for x in v_imp.p_imposable
            key = x[1];
            gen, ts, s = key;

            is_imp = false;
            if value(v_imp.p_is_imp[gen, ts, s]) > 0.5
                is_imp = true;
            end
            p_sol = value(v_imp.p_imposable[gen, ts, s]);
            p0 = launcher.uncertainties[gen, s, ts, ech];
            prev0 = launcher.previsions[gen, ts, ech];
            write(file, @sprintf("%s;%s;%s;%s;%s;%f;%f;%f\n", ech, gen, ts, s, is_imp, p_sol, p0, prev0));
        end
    end
end

function write_extra_output_csv(launcher::Workflow.Launcher, ech::DateTime, v_res::Workflow.ReserveModeler, v_flow::Dict{Tuple{String,DateTime,String},VariableRef})
    open(joinpath(launcher.dirpath, RESERVE_CSV), "a") do file
        if filesize(file) == 0
            write(file, @sprintf("%s;%s;%s;%s\n", "ech", "TS", "S", "res"));
        end

        dict_reserve_l = Dict{Tuple{DateTime,String}, Float64}()
        for ((ts_l,s_l), var_l) in v_res.p_res_pos
            if value(var_l) > 1e-6
                dict_reserve_l[ts_l, s_l] = value(var_l)
            end
        end
        for ((ts_l,s_l), var_l) in v_res.p_res_neg
            if value(var_l) > 1e-6
                @assert !haskey(dict_reserve_l, (ts_l,s_l)) #avoid having positive and negative reserves at once
                dict_reserve_l[ts_l, s_l] = -value(var_l) #negative reserve
            end
        end

        for ((ts_l,s_l), value_l) in dict_reserve_l
            write(file, @sprintf("%s;%s;%s;%f\n", ech, ts_l, s_l, value_l));
        end
    end

    open(joinpath(launcher.dirpath, FLOWS_CSV), "a") do file
        if filesize(file) == 0
            write(file, @sprintf("%s;%s;%s;%s;%s\n", "ech", "branch", "TS", "S", "flow"));
        end
        for ((branch_l, ts_l, s_l), flow_var_l) in v_flow
            write(file, @sprintf("%s;%s;%s;%s;%f\n", ech, branch_l, ts_l, s_l, value(flow_var_l)));
        end
    end
end

function write_previsions(launcher_p::Workflow.Launcher)
    open(joinpath(launcher_p.dirpath, SCHEDULE_CSV), "w") do file
        if filesize(file) == 0
            write(file, @sprintf("%s;%s;%s;%s;%s\n", "unit", "TS", "h_15m", "value", "decision_type"));
        end
        for ((gen_l, ts_l, ech_l), value_l) in launcher_p.previsions
            dmo_l = launcher_p.units[gen_l][5]
            decision_type_l = "flexible"
            if is_already_fixed(ech_l, ts_l, dmo_l)
                decision_type_l = "already_fixed"
            elseif is_to_decide(ech_l, ts_l, dmo_l)
                decision_type_l = "to_decide"
            end
            write(file, @sprintf("%s;%s;%s;%f;%s\n", gen_l, ts_l, ech_l, value_l, decision_type_l));
        end
    end
end
