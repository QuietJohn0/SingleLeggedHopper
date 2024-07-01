################## Double Pendulum 6.0 ##################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage as MF
using DifferentialEquations

##################### Start of Code #####################
### Create System ###
p = MF.System(0,1.6,.2,.135,.164,.013,.085,0,0,.085,1,40000,200,9.81,
            zeros(4, 2),0,zeros(2, 2),zeros(2, 2),zeros(2, 2),
            0,0,0,0,0,0,0,[0; 0],zeros(2, 2),
            0,0,0)


### Set up Fight Controller ###
p.setpoint = [0.6 -1.0; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### Set up Stance Controller ###
p.w = 50*pi
p.Fxd = 0
p.yhipd = ((p.Lᵦ + p.Lₐ)/3 + 2*(p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2]))/3)
p.Kp_s = [40; 20]
p.Kd_s = [30 10; 30 0]

### System Goals ###
p.dx_hop = .2
p.y_hop = .01/.05*p.dx_hop + .26
(p.Vxhipd,p.Vyhipd) = MF.VhipCalc(p)
(p.setpoint[2,1],p.setpoint[2,2]) = MF.Calc_Setpoint_Velcoity(p,.4)


### Initial Conditions ###
u_all = [0.60	-1.03	0.20	0.26	-3.75	4.83	-0.23	-1.20	0.00	-1.56]
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

(u_all,x_foot) = MF.invAug(p,u_all,10)

####### ODE Solver #######
tspan = (0, 5);

affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[3] + (u[6]<0)*-100 #-integrator.p.kₛ*u[3] - integrator.p.cₛ*u[6] + (u[6]<0)*100
cb_lift = ContinuousCallback(c_lift,affect!)

#u_all[end,:] = MF.Transition2Stanse(u_all[end,:],p)
(x₀,x_foot) = MF.invAug(p,u_all[end,:],6)
(p.A,p.ϕ) = MF.APhiGRFCalc(x₀,p)
prob = ODEProblem(MF.ode_S_fn!,x₀,tspan,p)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = .001)

##### Interpret Data #####
u_all = MF.q(u_all, p, sol.u[2:end,:], x_foot)
t_all = [[0]; sol.t[2:end,:]]

##### Export Data #####
name = "z_DataSet"
MF.ExportData(t_all,u_all,p,name)