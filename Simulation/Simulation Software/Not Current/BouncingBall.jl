#################### Double Pendulum ####################

################### Imported Librarys ###################
using DifferentialEquations
using Plots

####################### Functions #######################
### ODE Solver ###
function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -p[1]
  du[3] = u[4]
  du[4] = 0.0
end

### ODE Conditions ###
function condition(out,u,t,integrator) # Event when event_f(u,t) == 0
  out[1] = u[1]
  out[2] = (u[3] - 10.0)u[3]
end

### Condition Affects ###
function affect!(integrator, idx)
  display("boing")

  if idx == 1
    integrator.u[2] = -integrator.p[2]integrator.u[2]
  elseif idx == 2
    integrator.u[4] = -integrator.p[2]integrator.u[4]
  end

end

##################### Start of Code #####################
### Initial Conditions ###
u0 = [50.0,0.0,0.0,10.0]
tspan = (0.0,15.0)
p = [9.81, .7]

####### ODE Solver #######
prob = ODEProblem(f,u0,tspan,p)
cb = VectorContinuousCallback(condition,affect!,2)
sol = solve(prob,Tsit5(),callback=cb,dt=1e-1,adaptive = false)

####### Create Plot #######
plot(sol,idxs=(3,1))

####### Create GIF #######
n =  length(sol)
@gif for i ∈ 1:n
    scatter([sol.u[i][3]],[sol.u[i][1]],ylim = [0,50],xlim=[0,10],xlabel="X Position", ylabel="Y Position", 
        markercolors = [:red], shape = [:circle])
end every 1