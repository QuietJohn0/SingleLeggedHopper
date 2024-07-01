################## Double Pendulum Combined 9.0 ##################

################### Imported Librarys ###################
Path_name = "/Users/johna/OneDrive - Cal Poly/Documents/JULIACODE/MyRobotFunctionPackage/src"
if !(Path_name in LOAD_PATH)
    push!(LOAD_PATH, Path_name)
end
import .MyRobotFunctionPackage5 as MF5
using DifferentialEquations

##################### Start of Code #####################
### Create System ###
p = MF5.System(0.665,1.019,.128,.1362,.1636,.011,.0527,0.00207,0.00191,.00749,1,5000,300,40000,90,8.25,
            zeros(4, 2),zeros(2, 2),zeros(2, 2),
            0,0,0,0,0,0,zeros(2, 2), zeros(2, 2),zeros(2, 2),
            0,zeros(2, 2),
            0,0,0,0)

### System Goals ###
p.offset = .2
p.θ_liftOff = deg2rad(35)
p.Vxhipd = 1.25
p.Vyhipd = p.Vxhipd*tan(p.θ_liftOff)

### Set up Fight Controllers ###
p.setpoint = [0.72 -1.15; 0 0]
p.Kp_f = [35 0; 0 25]
p.Kd_f = [50 0; 0 10]

### Set up Stance Controller ###
p.w = 20*pi
p.Fxd = 0
p.Px = .6;
p.Py = .5;
p.Kp_s = [8 0; 0 8];
p.Kd_s = [0.04 0; 0 .04];
p.sp_s = MF5.InvKinSetpoint(p)

### Initial Conditions ###
u_all = [.6 -1.2 2.5-.06 .36 0 0 0 0 0 0]
(u_all,x_foot) = MF5.invAug(p,u_all,12)

t_cur = 0
t_all = [t_cur]

####### ODE Solver #######
tEnd = 10
savetime = .01

# Define Callbacks and Affect! Functions
affect!(integrator) = terminate!(integrator)
c_lift(u,t,integrator) = u[4] + (u[8]<0)*-100
c_impact(u,t,integrator) = MF5.getFootPos(u,integrator.p)[2] + (MF5.getFootVel(u,integrator.p)[2]>=0)*+100
cb_lift = ContinuousCallback(c_lift,affect!);
cb_impact = ContinuousCallback(c_impact,affect!);

###### Run Simulation ######
i = 1


while t_cur < tEnd -.01
    global t_cur, x₀, u_all, t_all, sw,i,x_foot,tspan,prob,sol, t_old
    tspan = (0, tEnd - t_cur)

    if i%2 ==1
        # Flight
        t_old = t_all[end]
        (p.setpoint[2,1],p.setpoint[2,2]) = MF5.Calc_Setpoint_Velcoity(p,u_all[end,:])
        (x₀,x_foot) = MF5.invAug(p,u_all[end,:],8)
        p.tf = MF5.FlightTimeApprox(x₀,p);
        p.C = MF5.TrajPlan(x₀,p.tf,p.setpoint);
        prob = ODEProblem(MF5.ode_F_fn!,x₀,tspan,p);
        sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_impact, abstol=1e-10,reltol=1e-10, saveat = savetime); 
    else
        # Stance
        (x₀,x_foot) = MF5.invAug(p,u_all[end,:],6)
        (p.A,p.ϕ) = MF5.APhiGRFCalc(x₀,p)
        prob = ODEProblem(MF5.ode_S_fn!,x₀,tspan,p)
        sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb_lift, abstol=1e-12,reltol=1e-12, saveat = savetime)
    end
    i += 1

    if length(sol.u[:,1]) >= 2 && sol.u[end-1,:] != sol.u[end,:]
        u_all = MF5.q(u_all, p, sol.u[2:end,:], x_foot)
        t_all = [t_all; ones(length(sol.t)-1,1)*t_cur + sol.t[2:end,:]]
        t_cur = t_all[end]
    elseif length(sol.u[:,1]) >= 3
        u_all = MF5.q(u_all, p, sol.u[2:end-1,:], x_foot)
        t_all = [t_all; ones(length(sol.t)-2,1)*t_cur + sol.t[2:end-1,:]]
        t_cur = t_all[end]
    else
        print("error")
    end  
end

##### Export Data #####
name = "z_DataSet"
MF5.ExportData(t_all,u_all,p,name)