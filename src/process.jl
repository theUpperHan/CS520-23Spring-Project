###############################
# Presolving
###############################

function remove_zero_columns(problem::IplpProblem)
    zero_columns = findall(sum(problem.A .!= 0, dims=1)[:] .== 0)

    # Deduce solution values x_i
    x_values = zeros(length(problem.c))
    for col in zero_columns
        if problem.c[col] < 0
            if problem.hi[col] == Inf
                return :unbounded
            else
                x_values[col] = problem.hi[col]
            end
        elseif problem.c[col] > 0
            x_values[col] = problem.lo[col]
        end
    end

    # Remove zero columns from A
    problem.A = problem.A[:, setdiff(1:size(problem.A, 2), zero_columns)]

    # Adjust the objective function coefficients and bounds
    problem.c = problem.c[setdiff(1:length(problem.c), zero_columns)]
    problem.lo = problem.lo[setdiff(1:length(problem.lo), zero_columns)]
    problem.hi = problem.hi[setdiff(1:length(problem.hi), zero_columns)]

    # Adjust the constraint vector b
    problem.b -= problem.A * x_values[setdiff(1:length(x_values), zero_columns)]

    return problem
end

function remove_zero_rows(problem::IplpProblem)
    zero_rows = findall(sum(problem.A .!= 0, dims=2)[:] .== 0)
    non_zero_rows = setdiff(1:size(problem.A, 1), zero_rows)
    problem.A = problem.A[non_zero_rows, :]
    problem.b = problem.b[non_zero_rows]
    return problem
end




###############################
# Standard Form
###############################

function convert_to_standard_form(P::IplpProblem)
    inf = 1.0e300
    m, n0 = size(P.A)

    # Initialize new matrices and vectors
    A_new = spzeros(m, 0)
    c_new = Vector{Float64}()

    for i in 1:n0
        lo_inf = P.lo[i] < -inf
        hi_inf = P.hi[i] > inf

        if lo_inf && hi_inf
            # x_i is unbounded: introduce x_i_pos and x_i_neg
            A_new = [A_new P.A[:, i] -P.A[:, i]]
            append!(c_new, [P.c[i], -P.c[i]])
        elseif lo_inf
            # x_i has an upper bound only: negate x_i
            A_new = [A_new -P.A[:, i]]
            append!(c_new, -P.c[i])
        else
            # x_i has a lower bound (or both): shift x_i
            A_new = [A_new P.A[:, i]]
            append!(c_new, P.c[i])
            P.b -= P.A[:, i] * P.lo[i]
        end
    end

    # Add slack variables for bounded variables
    idx_bounded = findall(x -> x > -inf && x < inf, P.hi)
    A_new = [A_new spzeros(m, length(idx_bounded))]
    c_new = vcat(c_new, zeros(length(idx_bounded)))
    I_bounded = Matrix{Float64}(I, length(idx_bounded), length(idx_bounded))

    As = [A_new; spzeros(length(idx_bounded), n0) I_bounded]
    bs = vcat(P.b, P.hi[idx_bounded] - P.lo[idx_bounded])
    cs = c_new

    lo_new = zeros(size(As, 2))
    hi_new = ones(size(As, 2)) * inf
    standard_form_ip = IplpProblem(cs, As, bs, lo_new, hi_new)

    return standard_form_ip
end


###############################
# Cholesky Factorization
###############################

function select_modified_cholesky_parameters(A::SparseMatrixCSC{Float64})
    beta = sqrt(cond(Matrix(A)))
    delta = eps(Float64) * norm(Matrix(A))
    return delta, beta
end

function modified_cholesky(A::SparseMatrixCSC{Float64}, delta::Real, beta::Real)
    @assert (n=size(A, 1)) == size(A, 2) "Input matrix must be a square matrix!"
    
    reorder = amd(A)
    L = A[reorder,reorder]

    d = ones(n)
    D = Diagonal(d)
        
    for i = 1:n-1
        theta = maximum(abs.(L[i+1:n, i]))
        d[i] = max(abs(L[i, i]), (theta / beta)^2, delta)
        L[i:n, i] ./= d[i]
        L[i, i] = 1.0
        L[i+1:n, i+1:n] .-= d[i] * (L[i+1:n, i] * L[i+1:n, i]')
    end

    d[n] = max(abs(L[n, n]), delta)
    L[n, n] = 1.0

    return (L=LowerTriangular(L), D=D, M=L, O=reorder)
end

function get_cholesky_lowtriangle(matrix, default=true)
    if default
        return cholesky(matrix)
    else
        delta, beta = select_modified_cholesky_parameters(matrix)
        F = modified_cholesky(matrix, delta, beta)
        return F.L
    end
end

