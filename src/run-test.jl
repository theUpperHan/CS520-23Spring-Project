include("test-helpers.jl")


if isempty(ARGS)
    println(ARGS)
    problems = ["lp_afiro","lp_brandy","lp_fit1d","lp_adlittle",
                "lp_agg","lp_ganges","lp_stocfor1", "lp_25fv47", "lpi_chemcom"]

    println("\nPerforming essential tests now....\n")
    for problem_name in problems
        try
            println("##################> Testing problem: ", problem_name, " <##################")
            D = test(problem_name, false)
            if D.op != -1
                println(GREEN, "Optimal Value: ", D.op, RESET)
            else
                println(RED, "No solution found,\nturn on info or see previous report", RESET)
            end
        catch e
            println(RED, "An error occurred: ", e, RESET)
        end
        println()
    end
    prompt_test()

else
    for filename in ARGS
        try
            D = test(filename)
            if D.op != -1
                println(GREEN, "Optimal Value: ", D.op, RESET)
            else
                println(RED, "No solution found,\nturn on info or see previous report", RESET)
            end
        catch e
            rintln(RED, "An error occurred: ", e, RESET)
        end
    end
end


# 