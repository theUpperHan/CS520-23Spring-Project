# Packages
using MatrixDepot
using SparseArrays
using AMD, LinearAlgebra
using Statistics 
using Printf
using HTTP
using Markdown

# Functions
include("prob-structs.jl")
include("process.jl")
include("solver-helpers.jl")
include("solver.jl")


function iplp(P::IplpProblem, info=true, tol=1e-8, maxit=100)
  P, zero_columns = remove_zero_columns(P)
  P, zero_rows = remove_zero_rows(P)
  As, bs, cs, bound1, bound2, bound3, bound4 = convert_to_standard_form(P)
  stdpb = IplpProblem(vec(cs), As, vec(bs), vec(P.lo), vec(P.hi))
  # Check Infeasibility
  @printf("===================Infeasible Check=====================\n")
  m,n = size(stdpb.A)
  A = [stdpb.A Matrix{Float64}(I,m,m)]
  b = copy(stdpb.b)
  c = [zeros(Float64, n); ones(Float64, m)]
  check_target = IplpProblem(c, A, b, vec(P.lo), vec(P.hi))
  check_x = predictor_corrector(check_target, false)[1]
  if abs(dot(c, check_x)) > 1e-8
      @printf("Infeasible LP Problem!\n")
      return IplpSolution(vec([0.]),false,vec(c),A,vec(b),vec([0.]),vec([0.]),vec([0.])), -1
  else
      @printf("Feasible LP Problem.\n")
  end
      
  if info
      @printf("\n=======================Solution=========================\n")
      @printf("The problem is feasible, finding solution now...\n")
      x_sol,λ_sol,s_sol,flag = predictor_corrector(stdpb)
  else
      x_sol,λ_sol,s_sol,flag = predictor_corrector(stdpb,false)
  end
  
  x = convert_x(P, zero_columns, bound1, bound2, bound3, bound4, x_sol)
  @show size(P.A)
  @show size(stdpb.A)
  @show size(x)
  @show size(x_sol)

  op_val = dot(P.c, x)
  if info
      @printf("Optimal Value: %.5f.\n", op_val)
  end
      
  return IplpSolution(vec(x),flag,vec(stdpb.c),stdpb.A,vec(stdpb.b),vec(x_sol),vec(λ_sol),vec(s_sol))
end

md = mdopen("LPnetlib/lpi_klein2")
pb = convert_matrixdepot(md)
sol = @time iplp(pb)