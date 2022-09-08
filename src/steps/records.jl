##########################################################################################
# Interface
##########################################################################################

abstract type RecordsWriter end

function get_destination(records_table::RecordsWriter)
    return records_table.destination
end

function set_destination!(records_table::RecordsWriter, destination_file_p::String)
    records_table.destination = destination_file_p
end

function write_header(records_table::RecordsWriter)
    output_file_l = get_destination(records_table)
    open(output_file_l, "w") do file_l
        Base.write(file_l, get_header(records_table))
        Base.write(file_l, "\n")
    end
end

function write_record(records_table::RecordsWriter)
    if get_config("EXTRA_LOG")
        println("WRITING :", formatted_record(records_table))
        output_file_l = get_destination(records_table)
        open(output_file_l, "a") do file_l
            Base.write(file_l, formatted_record(records_table))
            Base.write(file_l, "\n")
        end
    end
end

function write_record!(records_table::RecordsWriter)
    write_record(records_table)
    prepare_record!(records_table)
end

function formatted_record(records_table::RecordsWriter, sep=" ")
    return join(get_record(records_table), sep)
end

function get_header(records_table::RecordsWriter)::String
    error("interface needs to be implemented")
end

function get_record(records_table::RecordsWriter)::Vector{Any}
    error("interface needs to be implemented")
end

function init!(records_table::RecordsWriter, dir_p, filename_p, append=false)
    if get_config("EXTRA_LOG")
        mkpath(dir_p)
        set_destination!(records_table, joinpath(dir_p, filename_p))
        clean_record!(records_table)
        if (filesize(get_destination(records_table)) == 0) || !append
            write_header(records_table)
        end
    end
end

##############################################################################################
# Dynamic solve iterations Records
##############################################################################################
@with_kw mutable struct DynamicSolveRecords <: RecordsWriter
    destination::String = ""
    #Record: = 0
    iter::Int = 0
    nb_added_constraints::Int =-1
    nb_violated_constraints::Int =-1
    nb_total_constraints::Int =-1
    solver_time::Float64 = -1
    flow_computation_time::Float64 = -1
    flow_verification_time::Float64 = -1
    violations_filtering_time::Float64 = -1
end

function get_header(::DynamicSolveRecords)::String
    header_l = "iter nb_added_constraints nb_violated_constraints nb_total_constraints solver_time flow_computation_time flow_verification_time violations_filtering_time"
    return header_l
end

function get_record(records_table::DynamicSolveRecords)::Vector{Any}
    return [
        records_table.iter,
        records_table.nb_added_constraints,
        records_table.nb_violated_constraints,
        records_table.nb_total_constraints,
        records_table.solver_time,
        records_table.flow_computation_time,
        records_table.flow_verification_time,
        records_table.violations_filtering_time
    ]
end

function prepare_record!(records_table::DynamicSolveRecords)
    records_table.iter += 1
    records_table.nb_added_constraints = -1
    records_table.nb_violated_constraints = -1
    records_table.nb_total_constraints = -1
    records_table.solver_time = -1
    records_table.flow_computation_time = -1
    records_table.flow_verification_time = -1
    records_table.violations_filtering_time = -1
end

function clean_record!(records_table::DynamicSolveRecords)
    records_table.iter = 0
    records_table.nb_added_constraints = -1
    records_table.nb_violated_constraints = -1
    records_table.nb_total_constraints = -1
    records_table.solver_time = -1
    records_table.flow_computation_time = -1
    records_table.flow_verification_time = -1
    records_table.violations_filtering_time = -1
end

##############################################################################################
# TSO solve Records
##############################################################################################
@with_kw mutable struct TSOSolveRecords <: RecordsWriter
    destination::String = ""
    #Record: = 0
    usecase::String = ""
    dynamic::Bool = false
    nb_added_constraints::Int =-1
    nb_total_constraints::Int =-1
    nb_iters::Int = -1
    total_solving_time::Float64 = -1
    total_cstr_generation_time::Float64 = -1
end

function get_header(::TSOSolveRecords)::String
    header_l = "usecase dynamic nb_added_constraints nb_total_constraints nb_iters total_solving_time total_cstr_generation_time"
    return header_l
end

function get_record(records_table::TSOSolveRecords)::Vector{Any}
    return [
        records_table.usecase,
        records_table.dynamic,
        records_table.nb_added_constraints,
        records_table.nb_total_constraints,
        records_table.nb_iters,
        records_table.total_solving_time,
        records_table.total_cstr_generation_time
    ]
end

function prepare_record!(records_table::TSOSolveRecords)
    clean_record!(records_table)
end

function clean_record!(records_table::TSOSolveRecords)
    records_table.usecase = ""
    records_table.dynamic = false
    records_table.nb_added_constraints = -1
    records_table.nb_total_constraints = -1
    records_table.nb_iters = -1
    records_table.total_solving_time = -1
    records_table.total_cstr_generation_time = -1
end


##################
# Definitions
##################
const DYNAMIC_SOLVE_RECORDS = DynamicSolveRecords()
const TSO_SOLVE_RECORDS = TSOSolveRecords() #care if launching multiple TSOs in the same run, times of previous launches appear in the succeeding ones unless TIMER_TRACKS is cleared regularly
