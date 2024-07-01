################## Double Pendulum Combined 7.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage4 as MF4
using DifferentialEquations

##################### Start of Code #####################
### Create System ###
p = MF4.System(0,1.6,.2,.135,.164,.013,.085,0,0,.085,1,5000,250,40000,250,9.81,
            zeros(4, 2),zeros(2, 2),zeros(2, 2),
            0,0,0,0,0,0,zeros(2, 2), zeros(2, 2),zeros(2, 2),
            0,zeros(2, 2),
            0,0,0,0)

### Set up Fight Controllers ###
p.setpoint = [0.6 -1.0; -4 6]
p.Kp_f = [35 0; 0 25]
p.Kd_f = [50 0; 0 10]
### Set up Stance Controller ###
p.w = 50*pi
p.Fxd = 0
p.Px = .5;
p.Py = .5;
p.Kp_s = [.5 0; 0 .3];
p.Kd_s = [0.07 0; 0 .05];

### System Goals ###
p.offset = .5
p.θ_liftOff = deg2rad(60)
p.Vxhipd = 1.2
p.Vyhipd = p.Vxhipd*tan(p.θ_liftOff)

### Initial Conditions ###
u_all = [.6 -1.2 2.5-.06 .34 0 0 0 0 0 0]
(u_all,x_foot) = MF4.invAug(p,u_all,12)

t_cur = 0
t_all = [t_cur]

####### ODE Solver #######
tEnd = 6
savetime = .01

affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[4] + (u[8]<0)*-100
c_impact(u,t,integrator) = MF4.getFootPos(u,integrator.p)[2] + (MF4.getFootVel(u,integrator.p)[2]>=0)*+100
cb_lift = ContinuousCallback(c_lift,affect!);
cb_impact = ContinuousCallback(c_impact,affect!);

###### Run Simulation ######
i = 1


while t_cur < tEnd -.01
    global t_cur, x₀, u_all, t_all, sw,i,x_foot,tspan,prob,sol
    tspan = (0, tEnd - t_cur)

    if i%2 ==1
        (p.setpoint[2,1],p.setpoint[2,2]) = MF4.Calc_Setpoint_Velcoity(p,u_all[end,:])
        (x₀,x_foot) = MF4.invAug(p,u_all[end,:],8)
        p.tf = MF4.FlightTimeApprox(x₀,p);
        p.C = MF4.TrajPlan(x₀,p.tf,p.setpoint);
        prob = ODEProblem(MF4.ode_F_fn!,x₀,tspan,p);
        sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = savetime); 
    else
        (x₀,x_foot) = MF4.invAug(p,u_all[end,:],6)
        (p.A,p.ϕ) = MF4.APhiGRFCalc(x₀,p)
        (p.sp_s, p.tf) = MF4.InvKinSetpoint(p)
        p.C = MF4.TrajPlan(x₀,p.tf,p.sp_s);
        prob = ODEProblem(MF4.ode_S_fn!,x₀,tspan,p)
        sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = savetime)
    end
    i += 1

    if length(sol.u[:,1]) >= 2 && sol.u[end-1,:] != sol.u[end,:]
        u_all = MF4.q(u_all, p, sol.u[2:end,:], x_foot)
        t_all = [t_all; ones(length(sol.t)-1,1)*t_cur + sol.t[2:end,:]]
        t_cur = t_all[end]
    elseif length(sol.u[:,1]) >= 3
        u_all = MF4.q(u_all, p, sol.u[2:end-1,:], x_foot)
        t_all = [t_all; ones(length(sol.t)-2,1)*t_cur + sol.t[2:end-1,:]]
        t_cur = t_all[end]
    else
        print("bitty")
    end  
end

#=
for l in 1:length(t_all)
    println("time idx ",l,": ", t_all[l]," \tAM, H: ",MF4.Calculate_Angular_Momentum(p,u_all[l,:],true)[1]," \tU: ",u_all[l,5],", ",u_all[l,6])
end=#

##### Export Data #####
name = "z_DataSet"
MF4.ExportData(t_all,u_all,p,name)