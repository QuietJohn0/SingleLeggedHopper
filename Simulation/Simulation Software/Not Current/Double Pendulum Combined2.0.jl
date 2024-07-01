################## Double Pendulum Combined 1.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_fn!(du,u,p,t)
    if p[1] == 1
        s = p[2]
        τ =  [0 0 0 0]'
     
        M = [s.m₁*s.L₁^2 + s.m₂*(s.Lₐ^2 + s.L₂^2 + 2*s.Lₐ*s.L₂*cos(u[2]))      s.m₂*(s.L₂^2 + s.Lₐ*s.L₂*cos(u[2]))      (s.m₁*s.L₁+s.m₂*s.Lₐ)*cos(u[1]) + s.m₂*s.L₂*cos(u[1]+u[2])    (s.m₁*s.L₁+s.m₂*s.Lₐ)*sin(u[1]) + s.m₂*s.L₂*sin(u[1]+u[2]);
             s.m₂*(s.L₂^2 + s.Lₐ*s.L₂*cos(u[2]))                               s.m₂*s.L₂^2                               s.m₂*s.L₂*cos(u[1]+u[2])                                    s.m₂*s.L₂*sin(u[1]+u[2]);
            (s.m₁*s.L₁+s.m₂*s.Lₐ)*cos(u[1]) + s.m₂*s.L₂*cos(u[1]+u[2])          s.m₂*s.L₂*cos(u[1]+u[2])                   s.m₁ + s.m₂                                                0;
            (s.m₁*s.L₁+s.m₂*s.Lₐ)*sin(u[1]) + s.m₂*s.L₂*sin(u[1]+u[2])          s.m₂*s.L₂*sin(u[1]+u[2])                   0                                                          s.m₁ + s.m₂];
        V = [s.c₁*u[5] - s.m₂*s.Lₐ*s.L₂*u[6]*(u[6]+2*u[5])*sin(u[2]) + s.g*s.m₂*s.L₂*sin(u[1]+u[2]) + s.g*(s.m₁*s.L₁ + s.m₂*s.Lₐ)*sin(u[1]);
             s.c₂*u[6] + s.m₂*s.Lₐ*s.L₂*u[5]^2*sin(u[2])             + s.g*s.m₂*s.L₂*sin(u[1]+u[2]);
            -u[5]^2*(s.m₁*s.L₁ + s.m₂*s.Lₐ)*sin(u[2])      - (u[5]+u[6])^2*s.m₂*s.L₂*sin(u[1]+u[2]);
             u[5]^2*(s.m₁*s.L₁ + s.m₂*s.Lₐ)*cos(u[2])      + (u[5]+u[6])^2*s.m₂*s.L₂*cos(u[1]+u[2]) + s.g*(s.m₁ + s.m₂)];
        
        #fh = u[4] - s.Lₐ*cos(u[1]) - s.Lᵦ*cos(u[1]+u[2])

        q̈ = [M\(τ-V); 0]
        
    elseif p[1] == 2
        s = p[2]
        τ =  [0 0 0]'
        M = [s.Lᵦ^2*s.m₁ + 2*s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]) + (s.L₂-s.Lᵦ)^2*s.m₂ + (s.L₂-s.Lₐ)^2*s.m₁        s.Lᵦ^2*s.m₁ + (s.L₂-s.Lᵦ)^2*s.m₂ + s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2])        0;
             s.Lᵦ^2*s.m₁ +   s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]) + (s.L₂-s.Lᵦ)^2*s.m₂                             s.Lᵦ^2*s.m₁ + (s.L₂-s.Lᵦ)^2*s.m₂                                          0;
            -s.Lᵦ*s.m₁*sin(u[1] + u[2]) - (s.Lₐ-s.L₂)*s.m₁*sin(u[1]) - (s.Lᵦ - s.L₂)*s.m₂*sin(u[1] + u[2])    (-s.Lᵦ*s.m₁ - (s.Lᵦ - s.L₂)*s.m₂)*sin(u[1] + u[2])                        s.kₛ];

        V = [s.c₁*u[5] - u[6]*(2*u[5] + u[6])*s.Lᵦ*s.m₁*(s.Lₐ-s.L₂)*sin(u[2]) - s.m₁*s.g*(s.Lᵦ*sin(u[1] + u[2]) + (s.Lₐ-s.L₂)*sin(u[1])) - s.m₂*s.g*(s.Lᵦ-s.L₂)*sin(u[1] + u[2]);
             s.c₂*u[6] + u[5]^2*s.m₁*s.Lᵦ*(s.Lₐ-s.L₂)*sin(u[2])               - s.m₁*s.g*s.Lᵦ*sin(u[1] + u[2])                           - s.m₂*s.g*(s.Lᵦ-s.L₂)*sin(u[1] + u[2]);
            -((s.Lᵦ-s.L₂)*s.m₂ + s.Lᵦ*s.m₁)*(u[3]+u[4])^2*cos(u[1] + u[2]) - (s.Lₐ-s.L₂)*s.m₁*u[3]^2*cos(u[1]) + s.g*(s.m₁+s.m₂)]

        C = [1 0 0;
             0 1 0;
             -s.Lᵦ*cos(u[1]+u[2]) - s.Lₐ*cos(u[1])  -s.Lᵦ*cos(u[1]+u[2])    0;
             -s.Lᵦ*sin(u[1]+u[2]) - s.Lₐ*sin(u[1])  -s.Lᵦ*sin(u[1]+u[2])    1;
             0 0 1]
      
        D = [0
             0
             s.Lᵦ*(u[5]+u[6])^2*sin(u[1]+u[2]) + s.Lₐ*u[5]^2*sin(u[1]);
            -s.Lᵦ*(u[5]+u[6])^2*cos(u[1]+u[2]) - s.Lₐ*u[5]^2*cos(u[1]);
             0]
   
        q̈ = C*(M\(τ-V)) + D
        
    end

    du[1:4] = u[5:8]
    du[5:9] = q̈

    SA[du]
end

### ODE Conditions ###
function condition(out,u,t,integrator) # Event when event_f(u,t) == 0
    out[1] = u[4] - s.Lₐ*cos(u[1]) - s.Lᵦ*cos(u[1]+u[2])
    out[2] = (u[9] <  0) * u[9]
    if out[2] == 0
        #display(u[9])
    end
end

### Condition Affects ###
function affect!(integrator, idx)
    if idx == 1 && integrator.p[1] == 1
        display("hi")
        display(integrator.u[9])
        s = integrator.p[2]
        u = integrator.u
        A = [-s.Lᵦ*cos(u[1]+u[2]) - (s.Lₐ-s.L₁)*cos(u[1])   -s.Lᵦ*cos(u[1]+u[2]);
             -s.Lᵦ*sin(u[1]+u[2]) - (s.Lₐ-s.L₁)*sin(u[1])   -s.Lᵦ*sin(u[1]+u[2]);
             -(s.Lᵦ-s.L₂)*cos(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*cos(u[1]+u[2]) ;
             -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])]

        B = [s.L₁*cos(u[1])                       0                     1   0;
             s.L₁*sin(u[1])                       0                     0   1;
             s.Lₐ*cos(u[1]) + s.L₂*cos(u[1]+u[2]) s.L₂*cos(u[1]+u[2])   1   0;
             s.Lₐ*sin(u[1]) + s.L₂*sin(u[1]+u[2]) s.L₂*sin(u[1]+u[2])   0   1]

        M = [-m₁*(s.Lᵦ*cos(u[1]+u[2]) + (s.Lₐ-s.L₁)*cos(u[1]))  -m₁*(s.Lᵦ*sin(u[1]+u[2]) + (s.Lₐ-s.L₁)*sin(u[1]))   -m₂*((s.Lᵦ-s.L₂)*cos(u[1]+u[2]))    -m₂*((s.Lᵦ-s.L₂)*sin(u[1]+u[2]));
             -m₁*(s.Lₐ-s.L₁)*cos(u[1])                          -m₁*(s.Lₐ-s.L₁)*sin(u[1])                            0                                   0]

        C = [1 0;
             0 1;
            -s.Lᵦ*cos(u[1]+u[2]) - s.Lₐ*cos(u[1])  -s.Lᵦ*cos(u[1]+u[2]);
            -s.Lᵦ*sin(u[1]+u[2]) - s.Lₐ*sin(u[1])  -s.Lᵦ*sin(u[1]+u[2])]
        
        integrator.u[5:8] = C*inv(M*A)*M*B*u[5:8]
        integrator.p[1] = 2
    elseif idx == 2 && integrator.p[1] == 2
        display("hello")
        integrator.p[1] = 1
    end

end


### Convert to Cartesian ###
function X_Y(u,s)
     Lₐ,Lᵦ,L₁,L₂ = s
     xₒ = u[3,:]
     yₒ = u[4,:]
     u₁ = u[1,:]
     u₂ = u[2,:] + u[1,:]

     x = zeros(5,n)
     x[1,:] = xₒ
     x[2,:] = x[1,:] + L₁*sin.(u₁)
     x[3,:] = x[1,:] + Lₐ*sin.(u₁)
     x[4,:] = x[3,:] + L₂*sin.(u₂)
     x[5,:] = x[3,:] + Lᵦ*sin.(u₂)

     y = zeros(5,n)
     y[1,:] = yₒ
     y[2,:] = y[1,:] + -L₁*cos.(u₁)
     y[3,:] = y[1,:] + -Lₐ*cos.(u₁)
     y[4,:] = y[3,:] + -L₂*cos.(u₂)
     y[5,:] = y[3,:] + -Lᵦ*cos.(u₂)
     x,y
end

function PD_callback(integrator)
    u = integrator.u
    d = integrator.p[3].d
    if integrator.p[1] == 1
        τ = PDPosition(u,d,dt)
    elseif integrator.p[1] == 2
        τ = PDForceVcm(u,d,dt)
    end
    integrator.p[3].τk = τ
end

######## Controllers ########
function PDPosition(u,d)
    Kp = [2 4; 5 7]
    Kd = [2 4; 5 7]
    e = [u[1:4]-d[1:4]]
    τ = [Kp Kd]*e
end
function PDForceVcm(u,d,dt)
    Kp = [2 4; 5 7]
    Kd = [2 4; 5 7]
    e = [u[1:4]-d[1:4]]
    τ = [Kp Kd]*e
end

##### Convert Functions #####
function GRF2Torque(s,u,Bf)
    ## This function Will be used to convert desired GRF to expected torque
    M = [s.Lᵦ^2*s.m₁ + 2*s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]) + (s.L₂-s.Lᵦ)^2*s.m₂ + (s.L₂-s.Lₐ)^2*s.m₁        s.Lᵦ^2*s.m₁ + (s.L₂-s.Lᵦ)^2*s.m₂ + s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]);
         s.Lᵦ^2*s.m₁ +   s.Lᵦ*(s.Lₐ-s.L₂)*s.m₁*cos(u[2]) + (s.L₂-s.Lᵦ)^2*s.m₂                             s.Lᵦ^2*s.m₁ + (s.L₂-s.Lᵦ)^2*s.m₂];
    V = [s.c₁*u[5] - u[6]*(2*u[5] + u[6])*s.Lᵦ*s.m₁*(s.Lₐ-s.L₂)*sin(u[2]) - s.m₁*s.g*(s.Lᵦ*sin(u[1] + u[2]) + (s.Lₐ-s.L₂)*sin(u[1])) - s.m₂*s.g*(s.Lᵦ-s.L₂)*sin(u[1] + u[2]);
         s.c₂*u[6] + u[5]^2*s.m₁*s.Lᵦ*(s.Lₐ-s.L₂)*sin(u[2])               - s.m₁*s.g*s.Lᵦ*sin(u[1] + u[2])                           - s.m₂*s.g*(s.Lᵦ-s.L₂)*sin(u[1] + u[2])]
    A = [-s.Lᵦ*cos(u[1]+u[2]) - (s.Lₐ-s.L₁)*cos(u[1])   -s.Lᵦ*cos(u[1]+u[2]);
         -s.Lᵦ*sin(u[1]+u[2]) - (s.Lₐ-s.L₁)*sin(u[1])   -s.Lᵦ*sin(u[1]+u[2]);
         -(s.Lᵦ-s.L₂)*cos(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*cos(u[1]+u[2]) ;
         -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])                    -(s.Lᵦ-s.L₂)*sin(u[1]+u[2])]
    F = [ s.Lᵦ*(u[5]+u[6])^2*sin(u[1]+u[2]) + (s.Lₐ-s.L₁)*u[5]^2*sin(u[1]);
         -s.Lᵦ*(u[5]+u[6])^2*cos(u[1]+u[2]) - (s.Lₐ-s.L₁)*u[5]^2*cos(u[1]);
          (s.Lᵦ-s.L₂)*(u[5]+u[6])^2*sin(u[1]+u[2])
         -(s.Lᵦ-s.L₂)*(u[5]+u[6])^2*cos(u[1]+u[2])]
    mₘ = [s.m₁ 0    s.m₂ 0;
          0    s.m₁ 0    s.m₂]
    τ = M*inv(mₘ*A)*(Bf - mₘ*F - [0; g*(s.m₁ + s.m₂)]) + V    
end
function Vcm(s,u)
    J = [((s.m₁*s.L₁ + s.m₂*s.Lₐ)*cos(u[1]) + s.m₂*s.L₂*cos(u[1]+u[2]))/(s.m₁ + s.m₂)     s.m₂*s.L₂*cos(u[1]+u[2])/(s.m₁ + s.m₂)    1   0;
         ((s.m₁*s.L₁ + s.m₂*s.Lₐ)*sin(u[1]) + s.m₂*s.L₂*sin(u[1]+u[2]))/(s.m₁ + s.m₂)     s.m₂*s.L₂*sin(u[1]+u[2])/(s.m₁ + s.m₂)    0   1]
    Pcm = J*u[5:8]
    ẋcm = Pcm[1]
end

#function PI()

##################### Start of Code #####################
# State #
state = 1 # 1 Flight
          # 2 Stance
# System Parameters #
# [rad, rad,   m,  m, rad/s, rad/s, m/s, m/s] #
m₁,m₂,Lₐ,Lᵦ,L₁,L₂,c₁,c₂,kₛ,g = (4,2,.2,.18,.05,.06,.2,.2,1e8,9.81)
s = (m₁=m₁,m₂=m₂,Lₐ=Lₐ,Lᵦ=Lᵦ,L₁=L₁,L₂=L₂,c₁=c₁,c₂=c₂,kₛ=kₛ,g=g)
# Controller Stored Values #
tk,d,Vcmd,ωd,Ad,ϕd = (0, [-.4; -.2; 0; 0], .5, .3,50,.2)
c = (tk=tk,d=d,Vcmd=Vcmd,ωd=ωd,Ad=Ad,ϕd=ϕd)

p = [state s c]

### Initial Conditions ###
x₀ = [-.4, -.2, .24, .5, 0, 0, 0, -2, 0];

####### ODE Solver #######
tspan = (0.0, 7)
prob = ODEProblem(ode_fn!,x₀,tspan,p);
cb = VectorContinuousCallback(condition,affect!,2)
sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb,dt=1);

display("done")

##### Interpret Data #####
n =  length(sol)
dim = length(sol.u[1])
u = zeros(dim,n)
for i in 1:n
    u[:,i] = sol.u[i,1][1:dim]
end

w = Lₐ,Lᵦ,L₁,L₂
x,y = X_Y(u,w)

############### Phase Plot ###############
##### X V.S Y #####
#=plt = plot(x[1,:], y[1,:], lc=:black, lw=2)
plot!(x[3,:], y[3,:], lc=:red, lw=2)
plot!(x[5,:], y[5,:], lc=:blue, lw=2)
plot!([x[1,1]],[y[1,1]],markercolors = :red, shape = :circle, markersize = 3)
plot!([x[3,1]],[y[3,1]],markercolors = :red, shape = :circle, markersize = 3)
plot!([x[5,1]],[y[5,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Y and X Position of Orgin")
xlabel!("x")
ylabel!("y")
display(plt)=#

##### θ₁ V.S θ̇₁ #####
plt = plot(u[1,:], u[5,:], lc=:black, lw=2)
plot!([u[1,1]],[u[5,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₁ V.S θ̇₁")
xlabel!("θ₁")
ylabel!("θ̇₁")
display(plt)

##### θ₃ V.S θ̇₃ #####
plt = plot(u[2,:], u[6,:], lc=:black, lw=2)
plot!([u[2,1]],[u[6,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₃ V.S θ̇₃")
xlabel!("θ₁")
ylabel!("θ̇₃")
display(plt)

##### x V.S ẋ #####
plt = plot(u[3,:], u[7,:], lc=:black, lw=2)
plot!([u[3,1]],[u[7,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot x V.S ẋ")
xlabel!("x")
ylabel!("ẋ")
display(plt)

##### y V.S ẏ #####
plt = plot(u[4,:], u[8,:], lc=:black, lw=2)
plot!([u[4,1]],[u[8,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot y V.S ẏ")
xlabel!("y")
ylabel!("ẏ")
display(plt)


####### Create GIF #######
n =  length(sol)
s = p[2]
title = "Single Leg With State Transitions"
@gif for i ∈ 1:n
    #hi = sol.u[i][4] - s.Lₐ*cos(sol.u[i][1]) - s.Lᵦ*cos(sol.u[i][1]+ sol.u[i][2])
    hi = sol.u[i][4] - s.Lₐ*cos(sol.u[i][1]) - s.Lᵦ*cos(sol.u[i][1]+sol.u[i][2])
    plot(x[:,i],y[:,i],ylim = [-.5,.5],xlim=[-.8,.5],xlabel="X Position (m)", ylabel="Y Position (m)", 
       markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
       markersize = [3,8,6,8,16], linewidth = 4, title = hi-sol.u[i][9], aspect_ratio=:equal)
end every 1
