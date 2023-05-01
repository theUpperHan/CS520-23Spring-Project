include("src/test-helpers.jl")


if isempty(ARGS)
    problems = ["lp_afiro","lp_brandy","lp_fit1d","lp_adlittle",
                "lp_agg","lp_ganges","lp_stocfor1", "lp_25fv47", "lpi_chemcom"]

    problems_extra = ["lp_agg2", "lp_woodw", "lp_beaconfd", "lp_standata", "lp_degen2"
                        ,"lp_e226", "lp_grow15", "lpi_klein2", "lpi_vol1", "lpi_bgprtr"]

    println("\nPerforming essential tests now....\n")
    for problem_name in problems
        try
            println("##################> Testing problem: ", problem_name, " <##################")
            D = test(problem_name, false)
            if D.op != -1
                println(GREEN, "Optimal Value: ", D.op, RESET)
            else
                println(YELLOW, "No solution found due of infeasibility or needs more iterations\n", RESET)
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
            end
        catch e
            println(RED, "An error occurred: ", e, RESET)
        end
    end
end