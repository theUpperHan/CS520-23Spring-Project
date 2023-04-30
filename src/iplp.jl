# Packages
using MatrixDepot
using SparseArrays
using AMD, LinearAlgebra
using Statistics 
using Printf
using HTTP
using JuMP, GLPK

# Functions
include("prob-structs.jl")
include("process.jl")
include("solver-helpers.jl")
include("solver.jl")


function iplp(P::IplpProblem, info=true, tol=1e-8, maxit=100)
  P = remove_zero_columns(P)
  P = remove_zero_rows(P)
  stdpb = convert_to_standard_form(P)
  
  # Check Infeasibility
  if info
    @printf("===================Infeasible Check=====================\n")
  end

  m,n = size(stdpb.A)
  A = [stdpb.A Matrix{Float64}(I,m,m)]
  b = copy(stdpb.b)
  c = [zeros(Float64, n); ones(Float64, m)]
  check_target = IplpProblem(c, A, b, vec(P.lo), vec(P.hi))
  check_x = predictor_corrector(check_target, false)[1]
  if abs(dot(c, check_x)) > 1e-8
    if info
      @printf("Infeasible LP Problem!\n")
    end

    return IplpSolution(vec([0.]),false,vec(c),A,vec(b),vec([0.]),vec([0.]),vec([0.])), -1
  end
      
  if info
    @printf("Feasible LP Problem.\n")
    @printf("\n=======================Solution=========================\n")
    @printf("The problem is feasible, finding solution now...\n")
    x_sol,λ_sol,s_sol,flag = @time predictor_corrector(stdpb)
  else
    x_sol,λ_sol,s_sol,flag = @time predictor_corrector(stdpb,false)
  end
      
  op_val = dot(stdpb.c, x_sol)
  if info
    @printf("Optimal Value: %.5f.\n", op_val)
  end
      
  return IplpSolution(vec(x_sol),flag,vec(stdpb.c),stdpb.A,vec(stdpb.b),vec(x_sol),vec(λ_sol),vec(s_sol)), op_val
end

