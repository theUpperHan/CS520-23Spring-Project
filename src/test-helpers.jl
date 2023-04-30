include("iplp.jl")

function test(problem_name::String, info=true)
    md = mdopen("LPnetlib/" * problem_name)
    pb = convert_matrixdepot(md)

    solution, op_val = iplp(pb, info)
    
    if solution.flag
        @printf("Solution found.\n")
    else
        @printf("Not solution found after max iteration.\n")
    end
    
    if info
        @printf("\n======================Comparison========================\n")
    end
    
    A = md.A
    b = md.b
    c = md.c
    m, n = size(A)
    model = Model(with_optimizer(GLPK.Optimizer))
    
    @variable(model, x[1:n] >= 0)
    @objective(model, Min, sum(c[i] * x[i] for i in 1:n))
    for i in 1:m
        @constraint(model, sum(A[i, j] * x[j] for j in 1:n) <= b[i])
    end
    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        optimal_solution = value.(x)
        optimal_objective = objective_value(model)
        diff = optimal_objective-op_val
        absolute_diff = abs(diff)
        relative_diff = abs(diff/optimal_objective)
        if info
            println("GLPK Solution: ", optimal_objective)
            println("Absolute difference: ", absolute_diff)
            println("Relative difference: ", relative_diff)
        end
        return (abs=absolute_diff, rel=relative_diff, op=op_val)
    else
        println("No optimal solution found by GLPK solver.")
        return (abs=-1, rel=-1, op=op_val)
    end
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


