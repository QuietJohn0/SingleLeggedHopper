module MyRobotFunctionPackage4
########
# Implemented new control method for stance. using 
# inverse kinematics to convert desired hip position 2 angles.
########
################### Imported Librarys ###################
using LinearAlgebra
using Plots
using StaticArrays
using XLSX
using DataFrames

####################### Functions #######################
##### ODE Solver #####
    function ode_F_fn!(du,u,p,t)
        τ =  [Controller_F(u,p,t); 0; 0]

        m11 = p.I₁ + p.I₂ + p.L₁^2*p.m₁ + p.L₂^2*p.m₂ + 2*p.L₂*p.Lₐ*p.m₂*cos(u[2]) + p.Lₐ^2*p.m₂      
        m12 = p.I₂ + p.L₂^2*p.m₂ + p.L₂*p.Lₐ*p.m₂*cos(u[2])
        m13 = p.L₁*p.m₁*cos(u[1]) + p.L₂*p.m₂*cos(u[1] + u[2]) + p.Lₐ*p.m₂*cos(u[1])
        m14 = p.L₁*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) + p.Lₐ*p.m₂*sin(u[1])
        m21 = p.I₂ + p.L₂^2*p.m₂ + p.L₂*p.Lₐ*p.m₂*cos(u[2])
        m22 = p.I₂ + p.L₂^2*p.m₂
        m23 = p.L₂*p.m₂*cos(u[1] + u[2])
        m24 = p.L₂*p.m₂*sin(u[1] + u[2])
        m31 = p.L₁*p.m₁*cos(u[1]) + p.L₂*p.m₂*cos(u[1] + u[2]) + p.Lₐ*p.m₂*cos(u[1])
        m32 = p.L₂*p.m₂*cos(u[1] + u[2])
        m33 = p.m₀ + p.m₁ + p.m₂
        m34 = 0
        m41 = p.L₁*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) + p.Lₐ*p.m₂*sin(u[1])
        m42 = p.L₂*p.m₂*sin(u[1] + u[2])
        m43 = 0
        m44 = p.m₀ + p.m₁ + p.m₂        

        v1 = p.L₁*p.g*p.m₁*sin(u[1]) - 2*p.L₂*p.Lₐ*p.m₂*u[5]*u[6]*sin(u[2]) - p.L₂*p.Lₐ*p.m₂*u[6]^2*sin(u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.Lₐ*p.g*p.m₂*sin(u[1]) + p.c₁*u[5]
        v2 = p.L₂*p.Lₐ*p.m₂*u[5]^2*sin(u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.c₂*u[6]
        v3 = -p.L₁*p.m₁*u[5]^2*sin(u[1]) - p.L₂*p.m₂*u[5]^2*sin(u[1] + u[2]) - 2*p.L₂*p.m₂*u[5]*u[6]*sin(u[1] + u[2]) - p.L₂*p.m₂*u[6]^2*sin(u[1] + u[2]) - p.Lₐ*p.m₂*u[5]^2*sin(u[1])
        v4 = p.L₁*p.m₁*u[5]^2*cos(u[1]) + p.L₂*p.m₂*u[5]^2*cos(u[1] + u[2]) + 2*p.L₂*p.m₂*u[5]*u[6]*cos(u[1] + u[2]) + p.L₂*p.m₂*u[6]^2*cos(u[1] + u[2]) + p.Lₐ*p.m₂*u[5]^2*cos(u[1]) + p.g*p.m₀ + p.g*p.m₁ + p.g*p.m₂

        M = [m11 m12 m13 m14;
            m21 m22 m23 m24;
            m31 m32 m33 m34;
            m41 m42 m43 m44]
        
        V = [v1;
            v2;
            v3;
            v4]

        q̈ = M\(τ-V)

        du[1:4] = u[5:8]
        du[5:8] = q̈

        SA[du]
    end
    function ode_S_fn!(du,u,p,t)
        τ =  [Controller_S(u,p,t); 0; 0];

        v1 = 2*p.Lᵦ*p.L₁*p.m₁*u[5]*u[6]*sin(u[2]) + p.Lᵦ*p.L₁*p.m₁*u[6]^2*sin(u[2]) - 2*p.Lᵦ*p.Lₐ*p.m₁*u[5]*u[6]*sin(u[2]) - p.Lᵦ*p.Lₐ*p.m₁*u[6]^2*sin(u[2]) - p.Lᵦ*p.g*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.g*p.m₂*sin(u[1] + u[2]) + p.L₁*p.g*p.m₁*sin(u[1]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.g*p.m₁*sin(u[1]) + p.c₁*u[5]
        v2 = -p.Lᵦ*p.L₁*p.m₁*u[5]^2*sin(u[2]) + p.Lᵦ*p.Lₐ*p.m₁*u[5]^2*sin(u[2]) - p.Lᵦ*p.g*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.g*p.m₂*sin(u[1] + u[2]) + p.L₂*p.g*p.m₂*sin(u[1] + u[2]) + p.c₂*u[6]
        v3 =  p.Lᵦ*p.m₁*u[5]^2*sin(u[1] + u[2]) + 2*p.Lᵦ*p.m₁*u[5]*u[6]*sin(u[1] + u[2]) + p.Lᵦ*p.m₁*u[6]^2*sin(u[1] + u[2]) + p.Lᵦ*p.m₂*u[5]^2*sin(u[1] + u[2]) + 2*p.Lᵦ*p.m₂*u[5]*u[6]*sin(u[1] + u[2]) + p.Lᵦ*p.m₂*u[6]^2*sin(u[1] + u[2]) - p.L₁*p.m₁*u[5]^2*sin(u[1]) - p.L₂*p.m₂*u[5]^2*sin(u[1] + u[2]) - 2*p.L₂*p.m₂*u[5]*u[6]*sin(u[1] + u[2]) - p.L₂*p.m₂*u[6]^2*sin(u[1] + u[2]) + p.Lₐ*p.m₁*u[5]^2*sin(u[1]) + p.cₛ*u[7] + p.kₛₓ*u[3]
        v4 = -p.Lᵦ*p.m₁*u[5]^2*cos(u[1] + u[2]) - 2*p.Lᵦ*p.m₁*u[5]*u[6]*cos(u[1] + u[2]) - p.Lᵦ*p.m₁*u[6]^2*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*u[5]^2*cos(u[1] + u[2]) - 2*p.Lᵦ*p.m₂*u[5]*u[6]*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*u[6]^2*cos(u[1] + u[2]) + p.L₁*p.m₁*u[5]^2*cos(u[1]) + p.L₂*p.m₂*u[5]^2*cos(u[1] + u[2]) + 2*p.L₂*p.m₂*u[5]*u[6]*cos(u[1] + u[2]) + p.L₂*p.m₂*u[6]^2*cos(u[1] + u[2]) - p.Lₐ*p.m₁*u[5]^2*cos(u[1]) + p.cₛₓ*u[8] + p.g*p.m₁ + p.g*p.m₂ + p.kₛ*u[4]
        
        m11 = p.I₁ + p.I₂ + p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - 2*p.Lᵦ*p.L₁*p.m₁*cos(u[2]) - 2*p.Lᵦ*p.L₂*p.m₂ + 2*p.Lᵦ*p.Lₐ*p.m₁*cos(u[2]) + p.L₁^2*p.m₁ - 2*p.L₁*p.Lₐ*p.m₁ + p.L₂^2*p.m₂ + p.Lₐ^2*p.m₁
        m12 = p.I₂ + p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - p.Lᵦ*p.L₁*p.m₁*cos(u[2]) - 2*p.Lᵦ*p.L₂*p.m₂ + p.Lᵦ*p.Lₐ*p.m₁*cos(u[2]) + p.L₂^2*p.m₂
        m13 = -p.Lᵦ*p.m₁*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*cos(u[1] + u[2]) + p.L₁*p.m₁*cos(u[1]) + p.L₂*p.m₂*cos(u[1] + u[2]) - p.Lₐ*p.m₁*cos(u[1])
        m14 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₁*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.m₁*sin(u[1])
        m21 = p.I₂ + p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - p.Lᵦ*p.L₁*p.m₁*cos(u[2]) - 2*p.Lᵦ*p.L₂*p.m₂ + p.Lᵦ*p.Lₐ*p.m₁*cos(u[2]) + p.L₂^2*p.m₂
        m22 = p.I₂ + p.Lᵦ^2*p.m₁ + p.Lᵦ^2*p.m₂ - 2*p.Lᵦ*p.L₂*p.m₂ + p.L₂^2*p.m₂
        m23 = -p.Lᵦ*p.m₁*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*cos(u[1] + u[2]) + p.L₂*p.m₂*cos(u[1] + u[2])
        m24 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₂*p.m₂*sin(u[1] + u[2])
        m31 = -p.Lᵦ*p.m₁*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*cos(u[1] + u[2]) + p.L₁*p.m₁*cos(u[1]) + p.L₂*p.m₂*cos(u[1] + u[2]) - p.Lₐ*p.m₁*cos(u[1])
        m32 = -p.Lᵦ*p.m₁*cos(u[1] + u[2]) - p.Lᵦ*p.m₂*cos(u[1] + u[2]) + p.L₂*p.m₂*cos(u[1] + u[2])
        m33 = p.m₁ + p.m₂
        m34 = 0
        m41 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₁*p.m₁*sin(u[1]) + p.L₂*p.m₂*sin(u[1] + u[2]) - p.Lₐ*p.m₁*sin(u[1])
        m42 = -p.Lᵦ*p.m₁*sin(u[1] + u[2]) - p.Lᵦ*p.m₂*sin(u[1] + u[2]) + p.L₂*p.m₂*sin(u[1] + u[2])
        m43 = 0
        m44 = p.m₁ + p.m₂


        M = [m11 m12 m13 m14;
            m21 m22 m23 m24;
            m31 m32 m33 m34;
            m41 m42 m43 m44]
    
        V = [v1;
            v2;
            v3;
            v4]

        q̈ = M\(τ-V)

        du[1:4] = u[5:8]
        du[5:8] = q̈

        SA[du]
    end
### END ODE Solver ###

##### Controllers #####
    function Controller_F(u,p,t)
        θd = Traj(p,t)
        
        # Position θ₁
        θ₁d = θd[1,1]
        θ₁a = u[1]
        e_θ₁ = θ₁d - θ₁a

        # Position θ₃
        θ₃d = θd[1,2]
        θ₃a = u[2]
        e_θ₃ = θ₃d - θ₃a

        e = [e_θ₁;e_θ₃]

        # Position Derivative, θ̇₁
        θ̇₁d = θd[2,1]
        θ̇₁a = u[5]
        ė_θ₁ = θ̇₁d - θ̇₁a

        #  Position Derivative, θ̇₃
        θ̇₃d = θd[2,2]
        θ̇₃a = u[6]
        ė_θ₃ = θ̇₃d - θ̇₃a

        ė = [ė_θ₁;ė_θ₃]

        τ = p.Kp_f*e + p.Kd_f*ė
    end
    function Controller_S(u,p,t)
        # Feed Forward Torque
        τff = FeedForwardTorque(p,u,t)

        # Tragectorty Planning
        θd = Traj(p,t)
        
        # Position θ₁
        θ₁d = θd[1,1]
        θ₁a = u[1]
        e_θ₁ = θ₁d - θ₁a

        # Position θ₃
        θ₃d = θd[1,2]
        θ₃a = u[2]
        e_θ₃ = θ₃d - θ₃a

        e = [e_θ₁;e_θ₃]

        # Position Derivative, θ̇₁
        θ̇₁d = θd[2,1]
        θ̇₁a = u[5]
        ė_θ₁ = θ̇₁d - θ̇₁a

        #  Position Derivative, θ̇₃
        θ̇₃d = θd[2,2]
        θ̇₃a = u[6]
        ė_θ₃ = θ̇₃d - θ̇₃a

        ė = [ė_θ₁;ė_θ₃]

        τ = τff + p.Kp_s*e + p.Kd_s*ė
    end
### END Controllers ###

##### Flight Functions #####
    function FlightTimeApprox(u,p)
        yf = p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2])
        if((u[8]^2-2*p.g*(yf-u[4]))>=0)
            t = (u[8] + sqrt(u[8]^2-2*p.g*(yf-u[4])))/p.g
        else
            t = 0.005
        end
    end
    function Calc_Setpoint_Velcoity(p)
        dθ₁ = -2*p.Vxhipd*p.offset/(p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1] + p.setpoint[1,2]))
        dθ₃ = p.Vxhipd*p.offset/(p.Lᵦ*cos(p.setpoint[1,1] + p.setpoint[1,2]))
        (dθ₁,dθ₃)
    end
    function Calc_Setpoint_Velcoity(p,u)
        dθ₁ = -2*(p.Vxhipd-getFootVel(u,p)[1])*p.offset/(p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1] + p.setpoint[1,2]))
        dθ₃ = (p.Vxhipd-getFootVel(u,p)[1])*p.offset/(p.Lᵦ*cos(p.setpoint[1,1] + p.setpoint[1,2]))
        (dθ₁,dθ₃)
    end
### END Flight Functions ###

##### Stance Functions #####
    function Vcm(p,u)
    u = Aug(p,u,0)
    J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₀ + p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₀ + p.m₁ + p.m₂)    1   0;
            ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₀ + p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₀ + p.m₁ + p.m₂)    0   1]
    Pcm = J*u[5:8]
    ẋcm = Pcm[1]
    end
    function Acm(p,u)
    u = Aug(p,u,0)
    J = [((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2]))/(p.m₀ + p.m₁ + p.m₂)     p.m₂*p.L₂*cos(u[1]+u[2])/(p.m₀ + p.m₁ + p.m₂)    1   0;
            ((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]))/(p.m₀ + p.m₁ + p.m₂)     p.m₂*p.L₂*sin(u[1]+u[2])/(p.m₀ + p.m₁ + p.m₂)    0   1]
        
    Vⱼ= [-u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2])) + -u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*sin(u[1]+u[2]);
            u[5]^2*((p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])) +  u[6]*(u[6]+2*u[5])*p.m₂*p.L₂*cos(u[1]+u[2])]/(p.m₀ + p.m₁ + p.m₂)
    
    Pcm = J*u[5:8] + Vⱼ
    ẋcm = Pcm[1]
    end
    function Phip(p,u)
        xhip = Aug(p,u,0)[3:4]
    end
    function Vhip(p,u)
        ẋhip = Aug(p,u,0)[7:8]
    end
    function APhiGRFCalc(x₀,p)
    u = Aug(p,x₀,0)

    A = (((p.m₀+p.m₁+p.m₂)*p.w*(p.Vyhipd - u[8]))^2 + (p.cₛ*u[10])^2)/((p.m₀+p.m₁+p.m₂)*p.w*(p.Vyhipd - u[8]))
    ϕ = asin(-p.cₛ*u[10]/A)

    (A, ϕ)
    end
    function FeedForwardTorque(p,u,t)
        Fd = [p.Fxd; p.A*sin(p.w*t + p.ϕ)]
        J = [p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1])  p.Lᵦ*cos(u[1]+u[2]);
            p.Lᵦ*sin(u[1]+u[2]) + p.Lₐ*sin(u[1])  p.Lᵦ*sin(u[1]+u[2])]

        τ_mass = [p.m₁*p.g*p.L₁*sin(u[1]) + p.m₂*p.g*(p.L₂*sin(u[1]+u[2]) + p.Lₐ*sin(u[1]));
                  p.m₂*p.g*p.L₂*sin(u[1]+u[2])]
        τff = -J'*Fd + τ_mass
    end
    function InvKinSetpoint(p)  
        YD = (p.Lᵦ + p.Lₐ)*p.Py + (p.Lₐ*cos(p.setpoint[1,1]) + p.Lᵦ*cos(p.setpoint[1,1]+p.setpoint[1,2]))*(1-p.Py);
        XD = sqrt((p.Lᵦ + p.Lₐ)^2-YD^2)*p.Px;

        # Calculate Lift Off Joint Angles
        TH2_sign = sign(p.setpoint[1,2]);
        B = (YD^2 + XD^2 + p.Lₐ^2 - p.Lᵦ^2)/(-2*p.Lₐ);
        TH1 = 2*atan((XD+TH2_sign*sqrt(XD^2+YD^2-B^2))/(B-YD));
        TH2 = TH2_sign*acos((XD^2+YD^2-p.Lₐ^2-p.Lᵦ^2)/(2*p.Lᵦ*p.Lₐ));
    
        # Calculate Lift Off Joint Angular Velocities
        VXD = p.Vxhipd;
        VYD = p.Vyhipd;
        THD1 = (-VYD*p.Lᵦ*cos(TH1+TH2) + VXD*p.Lᵦ*sin(TH1+TH2))/(-p.Lᵦ*p.Lₐ*sin(TH2));
        THD2 = (VYD*(p.Lₐ*cos(TH1)+p.Lᵦ*cos(TH2+TH1)) - VXD*(p.Lₐ*sin(TH1)+p.Lᵦ*sin(TH1+TH2)))/(-p.Lᵦ*p.Lₐ*sin(TH2));
    
        # Define Setpoint
        SP = [TH1 TH2; THD1 THD2];
        TF = (pi-p.ϕ)/p.w;

        (SP, TF);
    end
### END Stance Functions ###

##### Trajectory Planning #####
    function TrajPlan(u,tf,setpoint)
        M_t =[1   0    0    0;
            0   1    0    0;
            1   tf   tf^2 tf^3;
            0   1    2*tf 3*tf^2]
        M_IC = [u[1:2]'; u[5:6]'; setpoint]
        C = (M_t^-1)*M_IC
    end
    function Traj(p,t)
        if p.tf < t
            t = p.tf
        end

        C = p.C
        θd = [1 t t^2 t^3;
            0 1 2*t 3*t^2]*C
    end
### END Trajectory Planning ###

##### Transition Functions #####
    function getFootPos(u,p)
    x_foot = [u[3] + p.Lₐ*sin(u[1]) + p.Lᵦ*sin(u[1]+ u[2]);
                u[4] - p.Lₐ*cos(u[1]) - p.Lᵦ*cos(u[1]+ u[2])]
    end
    function getFootVel(u,p)
    v_foot = [u[7] + u[5]*p.Lₐ*cos(u[1]) + (u[5]+ u[6])*p.Lᵦ*cos(u[1]+ u[2]);
                u[8] + u[5]*p.Lₐ*sin(u[1]) + (u[5]+ u[6])*p.Lᵦ*sin(u[1]+ u[2])]
    end
### END Transition Functions ###

##### Interpret Data #####
    function X_Y(u,p)
    n = length(u[:,1])

    xₒ = u[:,3]
    yₒ = u[:,4]
    u₁ = u[:,1]
    u₂ = u[:,2] + u[:,1]

    x = zeros(5,n)
    x[1,:] = xₒ
    x[2,:] = x[1,:] + p.L₁*sin.(u₁)
    x[3,:] = x[1,:] + p.Lₐ*sin.(u₁)
    x[4,:] = x[3,:] + p.L₂*sin.(u₂)
    x[5,:] = x[3,:] + p.Lᵦ*sin.(u₂)

    y = zeros(5,n)
    y[1,:] = yₒ
    y[2,:] = y[1,:] + -p.L₁*cos.(u₁)
    y[3,:] = y[1,:] + -p.Lₐ*cos.(u₁)
    y[4,:] = y[3,:] + -p.L₂*cos.(u₂)
    y[5,:] = y[3,:] + -p.Lᵦ*cos.(u₂)
    x,y
    end
    function Findsw(u_all)
        sw = [0]
        if (u_all[1,10] < 0)
            sw = [sw 1]
        end

        for i in 2:length(u_all[:,1])
            if (u_all[i,10] < 0 && u_all[i-1,10] >= 0)
                sw = [sw i]
            elseif (u_all[i,10] > 0 && u_all[i-1,10] <= 0)
                sw = [sw i-1]
            end
        end
        sw = [sw length(u_all[:,1])]
    end
    function findTorques(u_all,p,t_all, sw)
        τ = zeros(2,0)
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                #Stance
                (p.A,p.ϕ) = APhiGRFCalc(u_all[j,:],p)
                for l in j:k
                    τ = [τ Controller_S(u_all[l,:],p,t_all[l]-t_all[j])]
                end
            else
                #Flight
                (p.setpoint[2,1],p.setpoint[2,2]) = Calc_Setpoint_Velcoity(p,u_all[j,:])
                p.tf = FlightTimeApprox(u_all[j,:],p);
                for l in j:k
                    p.C = TrajPlan(u_all[j,:],p);
                    τ = [τ Controller_F(u_all[l,:],p,t_all[l]-t_all[j])]
                end
            end
        end
        τ
    end
    function Calculate_Angular_Momentum(p,u,isFlight)
        if (isFlight)
            ẋ = [u[5]; u[6]; u[7]; u[8]]
            A =[1                                  0                   0       0;
                1                                  1                   0       0;
                0                                  0                   1       0;
                0                                  0                   0       1;
                p.L₁*cos(u[1])                     0                   1       0;
                p.L₁*sin(u[1])                     0                   0       1;
                p.Lₐ*cos(u[1])+p.L₂*cos(u[1]+u[2]) p.L₂*cos(u[1]+u[2]) 1       0;
                p.Lₐ*sin(u[1])+p.L₂*sin(u[1]+u[2]) p.L₂*sin(u[1]+u[2]) 0       1]

        else
            ẋ = [u[5]; u[6]; u[11]; u[12]]
            A = [1                                               0                              0   0;
                1                                               1                              0   0;
                -p.Lᵦ*cos(u[1] + u[2]) - p.Lₐ*cos(u[1])        -p.Lᵦ*cos(u[1] + u[2])          1   0;
                -p.Lᵦ*sin(u[1] + u[2]) - p.Lₐ*sin(u[1])        -p.Lᵦ*sin(u[1] + u[2])          0   1;
                -p.Lᵦ*cos(u[1] + u[2]) - (p.Lₐ-p.L₁)*cos(u[1]) -p.Lᵦ*cos(u[1] + u[2])          1   0;
                -p.Lᵦ*sin(u[1] + u[2]) - (p.Lₐ-p.L₁)*sin(u[1]) -p.Lᵦ*sin(u[1] + u[2])          0   1;
                -(p.Lᵦ-p.L₂)*cos(u[1] + u[2])                  -(p.Lᵦ-p.L₂)*cos(u[1] + u[2])   1   0;
                -(p.Lᵦ-p.L₂)*sin(u[1] + u[2])                  -(p.Lᵦ-p.L₂)*sin(u[1] + u[2])   0   1]
        end

        m13 =       -p.m₀*p.Lₐ*cos(u[1])
        m14 =       -p.m₀*p.Lₐ*sin(u[1])
        m15 =       -p.m₁*(p.Lₐ-p.L₁)*cos(u[1])
        m16 =       -p.m₁*(p.Lₐ-p.L₁)*sin(u[1])
        m23 = m13 + -p.m₀*p.Lᵦ*cos(u[1]+u[2])
        m24 = m14 + -p.m₀*p.Lᵦ*sin(u[1]+u[2])
        m25 = m15 + -p.m₁*p.Lᵦ*cos(u[1]+u[2])
        m26 = m16 + -p.m₁*p.Lᵦ*sin(u[1]+u[2])    
        m27 =       -p.m₂*(p.Lᵦ-p.L₂)*cos(u[1]+u[2])
        m28 =       -p.m₂*(p.Lᵦ-p.L₂)*sin(u[1]+u[2])

        M = [p.I₁   0       m13 m14 m15 m16 0   0;
            p.I₁   p.I₂    m23 m24 m25 m26 m27 m28]
        
        H = M*A*ẋ
    end
    function Calculate_Angular_Momentum_All(p,u_all,sw,isFlight)
        H = zeros(2,0)
        for i in 1:length(sw)-1
            j = sw[i]+1
            k = sw[i+1]
            if i%2 == 0
                #Stance
                for l in j:k
                    H = [H Calculate_Angular_Momentum(p,u_all[l,:],isFlight)]
                end
            else
                #Flight
                for l in j:k
                    H = [H Calculate_Angular_Momentum(p,u_all[l,:],isFlight)]
                end
            end
        end
        H
    end
### END Interpret Data ###

##### State Mapping #####
    function q(u_all, p, u, x_foot)
        u_temp = reshape(Aug(p,u[1],x_foot),(1,12))
        for i in 2:length(u)
            u_temp = [u_temp; reshape(Aug(p,u[i],x_foot),(1,12))]
        end
        u_all = [u_all; u_temp]
    end
    function Aug(p,u,x_foot)
        if (length(u) == 8) && (u[4] <= 0)
            out = Aug_S(p,u,x_foot)
        elseif (length(u) == 8) && (u[4] > 0)
            out = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8] getFootPos(u,p)[1]-x_foot getFootPos(u,p)[2] getFootVel(u,p)[1] getFootVel(u,p)[2]]
        elseif length(u) == 12
            out = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8] u[9] u[10] u[11] u[12]]
        end
    end
    function invAug(p,u,n)
        if (n == 6)
            x₀ = [u[1] u[2] 0 getFootPos(u,p)[2] u[5] u[6] getFootVel(u,p)[1] getFootVel(u,p)[2]]
        elseif n == 8
            x₀ = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8]]
        elseif n == 12
            x₀ = [u[1] u[2] u[3] u[4] u[5] u[6] u[7] u[8] 0 getFootPos(u,p)[2] getFootVel(u,p)[1] getFootVel(u,p)[2]]
        end
        x_foot = getFootPos(u,p)[1]

        (x₀, x_foot)
    end
    function Aug_S(p,u,x_foot)
        C = [1 0 0 0;
             0 1 0 0;
            -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2])    1 0;
            -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])    0 1;]
        
        posȮ = C*u[5:8]
        
        posO = [x_foot + -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1]) + u[3]
                0      +  p.Lᵦ*cos(u[1]+u[2]) + p.Lₐ*cos(u[1]) + u[4]]

        u = [u[1] u[2] posO[1] posO[2] posȮ[1] posȮ[2] posȮ[3] posȮ[4] u[3] u[4] u[7] u[8]]
    end
### END State Mapping ###

##### Animating #####
    function PlotState!(x,y, title)
        plot!(x,y,ylim = [-.25,.85],xlim=[-1,1],
            xlabel="X Position (m)", 
            ylabel="Y Position (m)",
            title=title,
            markercolors = [:blue, :red, :blue, :red, :black], 
            markershape = [:circle,:circle, :circle, :circle, :circle],
            markersize = [3,5,5,5,3], 
            linewidth = 4, 
            aspect_ratio=:equal, 
            legend=false)
    end
    function PlotAnimation2(x_all,y_all,title,sw)
        x_all = AdjustXPositions(x_all)
        n =  length(x_all[1,:])
        anim = @animate for i ∈ 1:n
            plot([x_all[1,i]],[y_all[1,i]], markercolors=:black, markershapes=:rect, markersizes=:12)
            PlotState!(x_all[:,i],y_all[:,i],title)
            
            tot = length(sw)
            num = length(findall(sw.>i))
            if (tot - num) %2 == 0 
                plot!([x_all[5,i]-.15, x_all[5,i]+.15], [0,0], lc=:red, ls=:solid, lw=.7)
            else
                plot!([x_all[5,i]-.05, x_all[5,i]+.05], [0,0], lc=:green, ls=:solid, lw=.7)     
            end

        end every 1
        gif(anim, "anim_fps100.gif", fps = 100)
    end
    function PlotAnimation(x_all,y_all,title,sw)
        x_all = AdjustXPositions(x_all)
        n =  length(x_all[1,:])
        @gif for i ∈ 1:n
            plot([x_all[1,i]],[y_all[1,i]], markercolors=:black, markershapes=:rect, markersizes=:12)
            PlotState!(x_all[:,i],y_all[:,i],title)
            
            tot = length(sw)
            num = length(findall(sw.>i))
            if (tot - num) %2 == 0 
                plot!([x_all[5,i]-.15, x_all[5,i]+.15], [0,0], lc=:red, ls=:solid, lw=.7)
            else
                plot!([x_all[5,i]-.05, x_all[5,i]+.05], [0,0], lc=:green, ls=:solid, lw=.7)     
            end

        end every 1
    end
    function AdjustXPositions(x)
        for i in 1:length(x[1,:])
            if x[1,i] > 1.25 || x[1,i] < -1.25         
                n = (x[1,i] + 1.25) ÷ 2.5
                x[:,i] = x[:,i] .- n*2.5
            end
        end
        x
    end
### END Animating ###

##### Importing and Exporting Data #####
    function ExportData(t_all,u_all,p,name)
        df1 = DataFrame(t = t_all[:,1], θ₁=u_all[:,1], θ₃=u_all[:,2], x=u_all[:,3], y=u_all[:,4], θ̇₁=u_all[:,5], θ̇₃=u_all[:,6], ẋ=u_all[:,7], ẏ=u_all[:,8], xₛ=u_all[:,9], yₛ=u_all[:,10], ẋₛ=u_all[:,11], ẏₛ=u_all[:,12])
        df2 = DataFrame(setpoint1 = p.setpoint[:,1], setpoint2 = p.setpoint[:,2], Kp_f1 = p.Kp_f[:,1], Kp_f2 = p.Kp_f[:,2], Kd_f1 = p.Kd_f[:,1], Kd_f2 = p.Kd_f[:,2], 
                        w=[p.w, missing], Fxd=[p.Fxd, missing], Px=[p.Px, missing], Py=[p.Py, missing], Kp_s1 = p.Kp_s[:,1], Kp_s2 = p.Kp_s[:,2], Kd_s1 = p.Kd_s[:,1], Kd_s2 = p.Kd_s[:,2],
                        Vxhipd=[p.Vxhipd, missing], θ_liftOff=[p.θ_liftOff, missing],offset=[p.offset, missing])
        df3 = DataFrame(m₀=p.m₀,m₁=p.m₁,m₂=p.m₂,Lₐ=p.Lₐ,Lᵦ=p.Lᵦ,L₁=p.L₁,L₂=p.L₂,I₁=p.I₁,I₂=p.I₂,c₁=p.c₁,c₂=p.c₂,kₛₓ=p.kₛₓ,cₛₓ=p.cₛₓ,kₛ=p.kₛ,cₛ=p.cₛ,g=p.g)

        cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData")
        a = readdir()
        num = length(a)+1

        XLSX.writetable(name * string(num)*".xlsx", overwrite=true, 
            OUTPUT_DATA=(collect(DataFrames.eachcol(df1)), DataFrames.names(df1)),
            CONTROL_PARAM=(collect(DataFrames.eachcol(df2)), DataFrames.names(df2)),
            SYSTEM_PARAM=(collect(DataFrames.eachcol(df3)), DataFrames.names(df3)), 
        )
    end
    function ImportData()
        cd("C:\\Users\\johna\\OneDrive - Cal Poly\\Documents\\JavaVS-code\\CodeToTellAStory\\simulationData")
        a = readdir()
        filename = a[end]
        # Open the Excel file
        ImportData(filename)
    end
    function ImportData(filename)
        # Open the Excel file
        xf = XLSX.readxlsx(filename)
        sheetNames = XLSX.sheetnames(xf)

        # Extract from sheet 1
        sh = xf[sheetNames[1]]
        t_all = sh["A"][2:end]
        u_all = sh["B:M"][2:end,:]

        # Extract from sheet 2
        sh = xf[sheetNames[2]]
        p2 = sh["A:Q"][2:end,:]

        # Extract from sheet 3
        sh = xf[sheetNames[3]]
        p1 = sh["A:P"][2,:]

        # Build System
        p = System(p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8],p1[9],p1[10],p1[11],p1[12],p1[13],p1[14],p1[15],p1[16],
                [p2[:,1] p2[:,2]],[p2[:,3] p2[:,4]],[p2[:,5] p2[:,6]],
                0,p2[1,7],0,p2[1,8],p2[1,9],p2[1,10],zeros(2, 2), [p2[:,11] p2[:,12]], [p2[:,13] p2[:,14]],
                0, zeros(2, 2),
                p2[1,15],0,p2[1,16],p2[1,17])

        p.Vyhipd = p.Vxhipd*tan(p.θ_liftOff)

        (t_all,u_all,p,filename)
    end
### END Importing and Exporting Data ###

### System Structure ###
mutable struct System
    m₀::Float64
    m₁::Float64
    m₂::Float64
    Lₐ::Float64
    Lᵦ::Float64
    L₁::Float64
    L₂::Float64
    I₁::Float64
    I₂::Float64
    c₁::Float64
    c₂::Float64
    kₛₓ::Float64
    cₛₓ::Float64
    kₛ::Float64
    cₛ::Float64
    g::Float64

    setpoint::Matrix{Float64}
    Kp_f::Matrix{Float64}
    Kd_f::Matrix{Float64}

    A::Float64
    w::Float64
    ϕ::Float64
    Fxd::Float64
    Px::Float64
    Py::Float64
    sp_s::Matrix{Float64}
    Kp_s::Matrix{Float64}
    Kd_s::Matrix{Float64}

    tf::Float64
    C::Matrix{Float64}

    Vxhipd::Float64
    Vyhipd::Float64 
    θ_liftOff::Float64
    offset::Float64
end

hi() = print("sup")
greet() = print("Hello")

end # module MyRobotFunctionPackage3
