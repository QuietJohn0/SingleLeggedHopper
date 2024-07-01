################## Double Pendulum Combined 4.0 ##################

################### Imported Librarys ###################
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

### Set up Fight Controllers ###
p.setpoint = [0.8 -1.0; -4 6]
p.Kp_f = [25 15; 15 15]
p.Kd_f = [50 10; 10 10]

### Set up Stance Controller ###
p.w = 50*pi
p.Fxd = 0
p.yhipd = ((p.Lᵦ + p.Lₐ)/3 + 2*(p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2]))/3)
p.Kp_s = [20; 7]
p.Kd_s = [0 0; 10 0]
p.Kp_s = [0; 0]
p.Kd_s = [3 0; 10 0]

### System Goals ###
p.dx_hop = .2
p.y_hop = .01/.05*p.dx_hop + .26
p.offset = .4
(p.Vxhipd,p.Vyhipd) = MF.VhipCalc(p)

### Initial Conditions ###
u_all = [.6 -1.2 -.06 .34 0 0 0 0 0 0]
(u_all,x_foot) = MF.invAug(p,u_all,10)

t_cur = 0
t_all = [t_cur]

####### ODE Solver #######
tEnd = 3
savetime = .01

affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[3] + (u[6]<0)*-100
c_impact(u,t,integrator) = MF.getFootPos(u,integrator.p)[2] + (MF.getFootVel(u,integrator.p)[2]>=0)*+100
cb_lift = ContinuousCallback(c_lift,affect!)
cb_impact = ContinuousCallback(c_impact,affect!)

###### Run Simulation ######
i = 1
while t_cur < tEnd -.01
    global t_cur, x₀, u_all, t_all, sw,i
    tspan = (0, tEnd - t_cur)

    if i%2 ==1
        (p.setpoint[2,1],p.setpoint[2,2]) = MF.Calc_Setpoint_Velcoity(p)
        (x₀,x_foot) = MF.invAug(p,u_all[end,:],8)
        p.tf = MF.FlightTimeApprox(x₀,p);
        p.C = MF.TrajPlan(x₀,p);
        prob = ODEProblem(MF.ode_F_fn!,x₀,tspan,p);
        sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = savetime);    
    else
        u_all[end,:] = MF.Transition2Stanse(u_all[end,:],p)
        (x₀,x_foot) = MF.invAug(p,u_all[end,:],6)
        (p.A,p.ϕ) = MF.APhiGRFCalc(x₀,p)
        prob = ODEProblem(MF.ode_S_fn!,x₀,tspan,p)
        sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = savetime)
    end
    i += 1

    if length(sol.u[:,1]) >= 2 && sol.u[end-1,:] != sol.u[end,:]
        u_all = MF.q(u_all, p, sol.u[2:end,:], x_foot)
        t_all = [t_all; ones(length(sol.t)-1,1)*t_cur + sol.t[2:end,:]]
        t_cur = t_all[end]
    elseif length(sol.u[:,1]) >= 3
        u_all = MF.q(u_all, p, sol.u[2:end-1,:], x_foot)
        t_all = [t_all; ones(length(sol.t)-2,1)*t_cur + sol.t[2:end-1,:]]
        t_cur = t_all[end]
    else
        print("bitty")
    end
end

##### Export Data #####
name = "z_DataSet"
MF.ExportData(t_all,u_all,p,name)

