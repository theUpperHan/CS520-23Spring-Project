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

    return problem, zero_columns
end

function remove_zero_rows(problem::IplpProblem)
    zero_rows = findall(sum(problem.A .!= 0, dims=2)[:] .== 0)
    non_zero_rows = setdiff(1:size(problem.A, 1), zero_rows)
    problem.A = problem.A[non_zero_rows, :]
    problem.b = problem.b[non_zero_rows]
    return problem, zero_rows
end




###############################
# Standard Form
###############################

function convert_to_standard_form(P::IplpProblem)
    inf_val = 1.0e300

    rows, cols = size(P.A)

    unbounded_idx = Int64[]
    lower_bounded_idx = Int64[]
    upper_bounded_idx = Int64[]
    doubly_bounded_idx = Int64[]
    n = Int64[0, 0, 0, 0]

    for i = 1:cols
        if P.lo[i] < -inf_val
            if P.hi[i] > inf_val
                n[1] += 1
                push!(unbounded_idx, i)
            else
                n[3] += 1
                push!(upper_bounded_idx, i)
            end
        else
            if P.hi[i] > inf_val
                n[2] += 1
                push!(lower_bounded_idx, i)
            else
                n[4] += 1
                push!(doubly_bounded_idx, i)
            end
        end
    end

    c_standard = [P.c[unbounded_idx]; -P.c[unbounded_idx]; P.c[lower_bounded_idx];
                  -P.c[upper_bounded_idx]; P.c[doubly_bounded_idx]; zeros(n[4])]

    A_standard = [P.A[:, unbounded_idx] -P.A[:, unbounded_idx] P.A[:, lower_bounded_idx] -P.A[:, upper_bounded_idx] P.A[:, doubly_bounded_idx] spzeros(rows, n[4]);
                  spzeros(n[4], 2 * n[1] + n[2] + n[3]) Matrix{Float64}(I, n[4], n[4]) Matrix{Float64}(I, n[4], n[4])]

    b_standard = [P.b - P.A[:, lower_bounded_idx] * P.lo[lower_bounded_idx] - P.A[:, upper_bounded_idx] * P.hi[upper_bounded_idx] - P.A[:, doubly_bounded_idx] * P.lo[doubly_bounded_idx];
                  P.hi[doubly_bounded_idx] - P.lo[doubly_bounded_idx]]

    return A_standard, b_standard, c_standard, unbounded_idx, lower_bounded_idx, upper_bounded_idx, doubly_bounded_idx
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

