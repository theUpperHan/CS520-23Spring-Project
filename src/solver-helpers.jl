function starting_point(P::IplpProblem, default_cholesky=true)
    # 4 Steps
    # (1) Solving Primal and Dual starting points least squares
    # (2) Calculate delta primal and delta dual
    # (3) Calculate delta primal hat and delta dual hat
    # (4) Get starting point
    
    symm = P.A*P.A'
    if default_cholesky
        f = get_cholesky_lowtriangle(symm)
    else
        f = get_cholesky_lowtriangle(symm, false)
    end
    
    # (1)
    # Solve Primal
    d = f\P.b
    x_hat = P.A'*d
    
    # Solve Dual
    λ_hat = f\(P.A*P.c)
    s_hat = P.c - P.A'*λ_hat
    
    # (2)
    delta_p = max(-1.5*minimum(x_hat), 0)
    delta_d = max(-1.5*minimum(s_hat), 0)
    
    n = size(x_hat,1)
    e = ones(n,1)
    
    # (3)
    x_hat = x_hat + delta_p*e
    s_hat = s_hat + delta_d*e

    numerator = 0.5 * (x_hat'*s_hat)
    delta_p_hat = (numerator / (e'*s_hat))[1]
    delta_d_hat = (numerator / (e'*x_hat))[1]

    # (4)
    x_0 = x_hat + delta_p_hat*e
    λ_0 = λ_hat
    s_0 = s_hat + delta_d_hat*e

    return (x=x_0, λ=λ_0, s=s_0, f=f)
end

function build_kkt_matrix(A, x, s)
    m, n = size(A)
    zero_mn = spzeros(m, n)
    zero_nm = spzeros(n, m)
    zero_mm = spzeros(m, m)
    zero_nn = spzeros(n, n)

    X = spdiagm(0 => x[:, 1])
    Identity = Matrix{Float64}(I,n,n)
    S = spdiagm(0 => s[:, 1])

    M = [zero_mm  A         zero_mn;
         A'       zero_nn   Identity;
         zero_nm  S         X]

    return M
end

function solve_linear_system(A,x,s,rb,rc,rxs)
    KKT = build_kkt_matrix(A, x, s)
    f = lu(KKT)
    m = length(rb)
    n = length(rc)
    b = Array{Float64}([-rb; -rc; -rxs])
    b = f\b
    dlam = b[1:m]
    dx = b[1+m:m+n]
    ds = b[1+m+n:m+2*n]
    return dlam,dx,ds
end

function alpha_max(x,dx,hi=1.0)
    alpha = hi
    for i = 1:length(x)
        alpha = dx[i] < 0 ? min(alpha, -x[i] / dx[i]) : alpha
    end
    if alpha < 0
        alpha = Inf
    end
    alpha = min(alpha,hi)
    return alpha
end

function unbound_check(alpha_p, alpha_d)
    if alpha_p > 1e308 || alpha_d > 1e308
       return true 
    end
    return false
end

function select_gamma(P::IplpProblem)
    m, n = size(P.A)
    problem_size = max(m, n)
    
    if problem_size <= 100
        gamma = 0.01
    elseif problem_size <= 1000
        gamma = 0.05
    else
        gamma = 0.1
    end
    
    return gamma
end