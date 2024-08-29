################## Double Pendulum 7.0 ##################

################### Imported Librarys ###################
## Change the current working directory to the desired folder
Path_name = joinpath(@__DIR__, "MyRobotFunctionPackage/src")
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end

include(joinpath(Path_name, "MyRobotFunctionPackage5.jl"))
import .MyRobotFunctionPackage5 as MF5
using DifferentialEquations

##################### Start of Code #####################
### Create System ###
p = MF5.System(0.665,1.019,.128,.1362,.1636,.011,.0527,0.00207,0.00191,.00749,1,5000,300,40000,90,8.25,
            zeros(4, 2),zeros(2, 2),zeros(2, 2),
            0,0,0,0,0,0,zeros(2, 2), zeros(2, 2),zeros(2, 2),
            0,zeros(2, 2),
            0,0,0,0)

### Set up Fight Controllers ###
p.setpoint = [0.6 -1.0; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### System Goals ###
p.offset = .2
p.θ_liftOff = deg2rad(35)
p.Vxhipd = 1.25
p.Vyhipd = p.Vxhipd*tan(p.θ_liftOff)

### Initial Conditions ###
u_all = [.6 -1.2 .24 .34 0 0 .-.3 .3 0 0]
(u_all,x_foot) = MF5.invAug(p,u_all,12)

####### ODE Solver #######
tspan = (0,3);

affect!(integrator) = terminate!(integrator);
c_impact(u,t,integrator) = MF5.getFootPos(u,integrator.p)[2] + (MF5.getFootVel(u,integrator.p)[2]>=0)*+100
cb_impact = ContinuousCallback(c_impact,affect!);

(p.setpoint[2,1],p.setpoint[2,2]) = MF5.Calc_Setpoint_Velcoity(p,u_all[end,:])
(x₀,x_foot) = MF5.invAug(p,u_all[end,:],8)
p.tf = MF5.FlightTimeApprox(x₀,p);
p.C = MF5.TrajPlan(x₀,p.tf,p.setpoint);
prob = ODEProblem(MF5.ode_F_fn!,x₀,tspan,p);
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = .01);

##### Interpret Data #####
u_all = MF5.q(u_all, p, sol.u[2:end,:], x_foot)
t_all = [[0]; sol.t[2:end,:]]

##### Export Data #####
name = "z_DataSet"
MF5.ExportData(t_all,u_all,p,name)