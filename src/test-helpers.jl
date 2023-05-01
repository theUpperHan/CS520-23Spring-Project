include("iplp.jl")

function get_correct_ov(problem_name)
    if occursin("lpi", problem_name)
       return Inf
    end
    md_info = split(Markdown.plain(mdinfo("LPnetlib/" * problem_name)), "\n")
    header_line = "Name       Rows   Cols   Nonzeros    Bytes  BR      Optimal Value"
    header_idx = -1
    summary = ""
    for i in 1:length(md_info)
        if occursin(header_line, md_info[i])
            header_idx = i
            break
        end
    end
    
    summary = strip(md_info[header_idx+1], ['*',  ' '])
    pattern = r"(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S{0,1})\s+(-?\d+\.\d+E[+-]\d+)"
    m = match(pattern, summary)
    if m !== nothing
        parsed_dict = Dict(
            "Name" => m.captures[1],
            "Rows" => parse(Int, m.captures[2]),
            "Cols" => parse(Int, m.captures[3]),
            "Nonzeros" => parse(Int, m.captures[4]),
            "Bytes" => parse(Int, m.captures[5]),
            "BR" => m.captures[6],
            "OV" => parse(Float64, m.captures[7])
        )
    else
        println("No match found")
        parsed_dict = nothing
    end
    
    return parsed_dict["OV"]   
end

function test(problem_name::String, info=true)
    md = mdopen("LPnetlib/" * problem_name)
    pb = convert_matrixdepot(md)

    correct_ov = get_correct_ov(problem_name)

    res = @timed iplp(pb, info)
    
    solution, op_val = res.value
    elapsed_time = res.time
    println("IPLP solver takes: ", elapsed_time, " seconds to finish.")

    if solution.flag
        @printf("Solution found.\n")
    else
        @printf("Not solution found after max iteration.\n")
    end
    
    @printf("======================Comparison========================\n")
    
    if (correct_ov != Inf && !solution.flag) || (correct_ov == Inf && solution.flag)
        println(RED, "Not matched with correct result", RESET)
    else
        if correct_ov == Inf
            println(GREEN, "Successfully detect infeasibility", RESET)
        else
            println("Checking difference...")
            diff = correct_ov-op_val
            absolute_diff = abs(diff)
            relative_diff = abs(diff/correct_ov)
            println("Absolute difference: ", absolute_diff)
            println("Relative difference: ", relative_diff)
            
            if relative_diff < 1e-3
                println(GREEN, "Result is consistent with correct optimal value", RESET)
            else
                println(RED, "Result is wrong", GREEN)
            end
        end
    end

    return (cov=correct_ov, op=op_val)
end

function prompt_test(info=true)
    @printf("Which other test problems would you like to test? (Input as: lp_afiro, etc)\n")
    @printf("      [Type 'quit' or 'q' to quit]\n")
    problem_name = readline()

    if problem_name === nothing || problem_name == "" || problem_name == "quit" || problem_name == "q"
        return
    end

    test(problem_name, info)
end


