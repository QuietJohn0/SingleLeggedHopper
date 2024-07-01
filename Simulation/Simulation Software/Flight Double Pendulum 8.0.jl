################## Double Pendulum 7.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage3 as MF3
using DifferentialEquations

##################### Start of Code #####################
### Create System ###
p = MF3.System(0,1.6,.2,.135,.164,.013,.085,0,0,.085,1,5000,200,40000,200,9.81,
            zeros(4, 2),0,zeros(2, 2),zeros(2, 2),zeros(2, 2),
            0,0,0,0,0,0,0,[0; 0],zeros(2, 2),
            0,0,0)

### Set up Fight Controllers ###
p.setpoint = [0.6 -1.0; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### System Goals ###
p.dx_hop = .2
p.y_hop = .01/.05*p.dx_hop + .26
p.offset = .4
(p.Vxhipd,p.Vyhipd) = MF3.VhipCalc(p)

### Initial Conditions ###
u_all = [.6 -1.2 .24 .34 0 0 .-.3 .3 0 0]
(u_all,x_foot) = MF3.invAug(p,u_all,12)

####### ODE Solver #######
tspan = (0,3);

affect!(integrator) = terminate!(integrator);
c_impact(u,t,integrator) = MF3.getFootPos(u,integrator.p)[2] + (MF3.getFootVel(u,integrator.p)[2]>=0)*+100
cb_impact = ContinuousCallback(c_impact,affect!);

(p.setpoint[2,1],p.setpoint[2,2]) = MF3.Calc_Setpoint_Velcoity(p,u_all[end,:])
(x₀,x_foot) = MF3.invAug(p,u_all[end,:],8)
p.tf = MF3.FlightTimeApprox(x₀,p);
p.C = MF3.TrajPlan(x₀,p);
prob = ODEProblem(MF3.ode_F_fn!,x₀,tspan,p);
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = .01);

##### Interpret Data #####
u_all = MF3.q(u_all, p, sol.u[2:end,:], x_foot)
t_all = [[0]; sol.t[2:end,:]]

##### Export Data #####
name = "z_DataSet"
MF3.ExportData(t_all,u_all,p,name)