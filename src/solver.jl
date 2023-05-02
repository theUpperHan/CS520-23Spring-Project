function predictor_corrector(P::IplpProblem, info=true, tol=1e-8, maxit=100)    
    iter = 0
    SP = starting_point(P)
    GAMMA = select_gamma(P) #limit step length
    m,n = size(P.A)
    
    x_0 = copy(SP.x)
    λ_0 = copy(SP.λ)
    s_0 = copy(SP.s)
    
    if info
        @printf("%10s %12s %12s\n", "Iteration", "Alpha Prim", "Alpha Dual")
        @printf("%10d %12.4g %12.4g\n", iter, 0., 0.)
    end
    
    x1,λ1,s1,flag = vec([0.]),vec([0.]),vec([0.]),false
    
    for iter=1:maxit
        # Residuals for primal feasibility       rb
        #               dual feasibility         rc
        #               complementary slackness  rxs
        rb  = P.A*x_0-P.b
        rc  = P.A'*λ_0+s_0-P.c
        rxs = x_0.*s_0

        # Predictor
        λ_aff,x_aff,s_aff = solve_linear_system(P.A,x_0,s_0,rb,rc,rxs) # direction
        alpha_aff_primal  = alpha_max(x_0,x_aff) # step lengths
        alpha_aff_dual = alpha_max(s_0,s_aff)
        
        # Corrector
        μ = mean(rxs)
        μ_affine = sum((x_0 .+ alpha_aff_primal .* x_aff) .* (s_0 .+ alpha_aff_dual .* s_aff)) / n
        σ = (μ_affine/μ)^3 # centering parameter
        rxs = x_aff.*s_aff.-σ*μ
        λ_cc,x_cc,s_cc = solve_linear_system(P.A,x_0,s_0,spzeros(m),spzeros(n),rxs)

        # Combine to find search direction
        direction_x = x_aff+x_cc
        direction_λ = λ_aff+λ_cc
        direction_s = s_aff+s_cc
        
        # Finding the max step length it can take 
        # Find max allowable step len Wright Chapter 4
        alpha_max_pri, alpha_max_dual = alpha_max(x_0, direction_x, Inf), alpha_max(s_0, direction_s, Inf)
        
        # Wright P205 Chapter 10
        x1_pri, s1_dual = x_0 + alpha_max_pri * direction_x, s_0 + alpha_max_dual * direction_s
        μp= dot(x1_pri, s1_dual) / n
        f_prim = (GAMMA * μp / s1_dual[argmin(x1_pri)] - x_0[argmin(x1_pri)]) / alpha_max_pri / direction_x[argmin(x1_pri)]
        f_dual = (GAMMA * μp / x1_pri[argmin(s1_dual)] - s_0[argmin(s1_dual)]) / alpha_max_dual / direction_s[argmin(s1_dual)]
        alpha_pri, alpha_dual = max(1 - GAMMA, f_prim) * alpha_max_pri, max(1 - GAMMA, f_dual) * alpha_max_dual

        
        # Unboundness check
        if unbound_check(alpha_pri, alpha_dual)
            warn("This problem is unbounded")
            return x1,λ1,s1,false,iter
        end
        
        # Take the step
        x1 = x_0+alpha_pri*direction_x
        λ1 = λ_0+alpha_dual*direction_λ
        s1 = s_0+alpha_dual*direction_s
        x_0 = x1
        λ_0 = λ1
        s_0 = s1
        
        if info
            @printf("%10d %12.4g %12.4g\n", iter, alpha_pri, alpha_dual)
        end
        
        r1 = norm(P.A*x1-P.b)/(1+norm(P.b))
        r2 = norm(P.A'*λ1+s1-P.c)/(1+norm(P.c))
        r3 = abs(dot(P.c,x1)-dot(P.b,λ1))/(1+abs(dot(P.c,x1)))
        if r1 < tol && r2 < tol && r3 < tol
            flag = true
            break
        end
    end
    
    return x1,λ1,s1,flag
end