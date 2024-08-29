################## Double Pendulum 6.0 ##################
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


### Set up Fight Controller ###
p.setpoint = [0.6 -1.0; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### Set up Stance Controller ###
p.w = 20*pi
p.Fxd = 0
p.Px = .6;
p.Py = .5;
p.Kp_s = [8 0; 0 8];
p.Kd_s = [0.04 0; 0 .04];
p.sp_s = MF5.InvKinSetpoint(p)

### System Goals ###
p.offset = .2
p.θ_liftOff = deg2rad(35)
p.Vxhipd = 1.25
p.Vyhipd = p.Vxhipd*tan(p.θ_liftOff)


### Initial Conditions ###
u_all = [0.8089202258561855
-1.0238328067600322
-0.06753826499825777
 0.2534149998937212
-7.160650702004575
 5.272175239167533
 0.14868554594874597
-1.1532823714820515
 5.551115123125783e-17
-1.786672494037899]

(u_all,x_foot) = MF5.invAug(p,u_all,12)

####### ODE Solver #######
tspan = (0, 5);

affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[4] + (u[8]<0)*-100
cb_lift = ContinuousCallback(c_lift,affect!)

(x₀,x_foot) = MF5.invAug(p,u_all[end,:],6)
(p.A,p.ϕ) = MF5.APhiGRFCalc(x₀,p)
prob = ODEProblem(MF5.ode_S_fn!,x₀,tspan,p)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = .001)

##### Interpret Data #####
u_all = MF5.q(u_all, p, sol.u[2:end,:], x_foot)
t_all = [[0]; sol.t[2:end,:]]

##### Export Data #####
name = "z_DataSet"
MF5.ExportData(t_all,u_all,p,name)