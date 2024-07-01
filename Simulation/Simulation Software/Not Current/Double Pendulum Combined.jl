################## Double Pendulum Combined 1.0 ##################

################### Imported Librarys ###################
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

####################### Functions #######################
### ODE Solver ###
function ode_fn!(du,u,p,t)
    if p[2] == 1
        p = p[1]
        τ =  [0 0 0 0]'
     
        M = [p.m₁*p.L₁^2 + p.m₂*(p.Lₐ^2 + p.L₂^2 + 2*p.Lₐ*p.L₂*cos(u[2]))      p.m₂*(p.L₂^2 + p.Lₐ*p.L₂*cos(u[2]))      (p.m₁*p.L₁+p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])    (p.m₁*p.L₁+p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2]);
            p.m₂*(p.L₂^2 + p.Lₐ*p.L₂*cos(u[2]))                               p.m₂*p.L₂^2                               p.m₂*p.L₂*cos(u[1]+u[2])                                    p.m₂*p.L₂*sin(u[1]+u[2]);
            (p.m₁*p.L₁+p.m₂*p.Lₐ)*cos(u[1]) + p.m₂*p.L₂*cos(u[1]+u[2])          p.m₂*p.L₂*cos(u[1]+u[2])                   p.m₁ + p.m₂                                                0;
            (p.m₁*p.L₁+p.m₂*p.Lₐ)*sin(u[1]) + p.m₂*p.L₂*sin(u[1]+u[2])          p.m₂*p.L₂*sin(u[1]+u[2])                   0                                                          p.m₁ + p.m₂];
        V = [p.c₁*u[5] - p.m₂*p.Lₐ*p.L₂*u[6]*(u[6]+2*u[5])*sin(u[2]) + p.g*p.m₂*p.L₂*sin(u[1]+u[2]) + p.g*(p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[1]);
            p.c₂*u[6] + p.m₂*p.Lₐ*p.L₂*u[5]^2*sin(u[2])             + p.g*p.m₂*p.L₂*sin(u[1]+u[2]);
            -u[5]^2*(p.m₁*p.L₁ + p.m₂*p.Lₐ)*sin(u[2])      - (u[5]+u[6])^2*p.m₂*p.L₂*sin(u[1]+u[2]);
            u[5]^2*(p.m₁*p.L₁ + p.m₂*p.Lₐ)*cos(u[2])      + (u[5]+u[6])^2*p.m₂*p.L₂*cos(u[1]+u[2]) + p.g*(p.m₁ + p.m₂)];
        
        q̈ = M\(τ-V)
        
    elseif p[2] == 2
        p = p[1]
        τ =  [0 0]'
        M = [p.Lᵦ^2*p.m₁ + 2*p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]) + (p.L₂-p.Lᵦ)^2*p.m₂ + (p.L₂-p.Lₐ)^2*p.m₁        p.Lᵦ^2*p.m₁ + (p.L₂-p.Lᵦ)^2*p.m₂ + p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]);
             p.Lᵦ^2*p.m₁ +   p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]) + (p.L₂-p.Lᵦ)^2*p.m₂                             p.Lᵦ^2*p.m₁ + (p.L₂-p.Lᵦ)^2*p.m₂];
       
        V = [p.c₁*u[5] - u[6]*(2*u[5] + u[6])*p.Lᵦ*p.m₁*(p.Lₐ-p.L₂)*sin(u[2]) - p.m₁*p.g*(p.Lᵦ*sin(u[1] + u[2]) + (p.Lₐ-p.L₂)*sin(u[1])) - p.m₂*p.g*(p.Lᵦ-p.L₂)*sin(u[1] + u[2]);
             p.c₂*u[6] + u[5]^2*p.m₁*p.Lᵦ*(p.Lₐ-p.L₂)*sin(u[2])               - p.m₁*p.g*p.Lᵦ*sin(u[1] + u[2])                           - p.m₂*p.g*(p.Lᵦ-p.L₂)*sin(u[1] + u[2])]
                                                                               
        C = [1 0;
             0 1;
             -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2]);
             -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])]
      
        D = [0
             0
             p.Lᵦ*(u[5]+u[6])^2*sin(u[1]+u[2]) + p.Lₐ*u[5]^2*sin(u[1]);
            -p.Lᵦ*(u[5]+u[6])^2*cos(u[1]+u[2]) - p.Lₐ*u[5]^2*cos(u[1])]
   
        q̈ = C*(M\(τ-V)) + D
        
    end

    du[1:4] = u[5:8]
    du[5:8] = q̈

    SA[du]
end

### ODE Conditions ###
function condition(out,u,t,integrator) # Event when event_f(u,t) == 0
    p = integrator.p[1]
    τ = [0 0]'
    out[1] = 1
    out[2] = 1

    M = [p.Lᵦ^2*p.m₁ + 2*p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]) + (p.L₂-p.Lᵦ)^2*p.m₂ + (p.L₂-p.Lₐ)^2*p.m₁        p.Lᵦ^2*p.m₁ + (p.L₂-p.Lᵦ)^2*p.m₂ + p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]);
         p.Lᵦ^2*p.m₁ +   p.Lᵦ*(p.Lₐ-p.L₂)*p.m₁*cos(u[2]) + (p.L₂-p.Lᵦ)^2*p.m₂                             p.Lᵦ^2*p.m₁ + (p.L₂-p.Lᵦ)^2*p.m₂];
    V = [p.c₁*u[5] - u[6]*(2*u[5] + u[6])*p.Lᵦ*p.m₁*(p.Lₐ-p.L₂)*sin(u[2]) - p.m₁*p.g*(p.Lᵦ*sin(u[1] + u[2]) + (p.Lₐ-p.L₂)*sin(u[1])) - p.m₂*p.g*(p.Lᵦ-p.L₂)*sin(u[1] + u[2]);
         p.c₂*u[6] + u[5]^2*p.m₁*p.Lᵦ*(p.Lₐ-p.L₂)*sin(u[2])               - p.m₁*p.g*p.Lᵦ*sin(u[1] + u[2])                           - p.m₂*p.g*(p.Lᵦ-p.L₂)*sin(u[1] + u[2])]
    A = [-p.Lᵦ*sin(u[1]+u[2]) - (p.Lₐ-p.L₁)*sin(u[1])   -p.Lᵦ*sin(u[1]+u[2]);
         -(p.Lᵦ-p.L₂)*sin(u[1]+u[2])                    -(p.Lᵦ-p.L₂)*sin(u[1]+u[2])]
    F = [-p.Lᵦ*(u[5]+u[6])^2*cos(u[1]+u[2]) - (p.Lₐ-p.L₁)*u[5]^2*cos(u[1]);
        -(p.Lᵦ-p.L₂)*(u[5]+u[6])^2*cos(u[1]+u[2])]
    GRF = ([p.m₁ p.m₂]*(p.g .+ A*(M\(τ-V)) + F))[1]
    pos = (u[4] - p.Lₐ*cos(u[1]) - p.Lᵦ*cos(u[1]+u[2]))
    if (integrator.p[2] == 1) & (GRF > 0)
        out[1] = pos
        display(string(pos) *", "* string(GRF))
    end
    if (integrator.p[2] == 2)
        out[2] = GRF 
        #display(string(GRF) *", " * string(pos))
        # out[2] = 1.0 * !((p.g*(p.m₁+p.m₂) + ([0 0 0 p.m₁ 0 p.m₂]*(A*u[1:2] + F))[1] < 0) && (integrator.p[2] == 2))
    end
end

### Condition Affects ###
function affect!(integrator, idx)
    if idx == 1
        display("hi")
        p = integrator.p[1]
        u = integrator.u
        A = [-p.Lᵦ*cos(u[1]+u[2]) - (p.Lₐ-p.L₁)*cos(u[1])   -p.Lᵦ*cos(u[1]+u[2]);
             -p.Lᵦ*sin(u[1]+u[2]) - (p.Lₐ-p.L₁)*sin(u[1])   -p.Lᵦ*sin(u[1]+u[2]);
             -(p.Lᵦ-p.L₂)*cos(u[1]+u[2])                    -(p.Lᵦ-p.L₂)*cos(u[1]+u[2]) ;
             -(p.Lᵦ-p.L₂)*sin(u[1]+u[2])                    -(p.Lᵦ-p.L₂)*sin(u[1]+u[2])]

        B = [p.L₁*cos(u[1])                       0                     1   0;
             p.L₁*sin(u[1])                       0                     0   1;
             p.Lₐ*cos(u[1]) + p.L₂*cos(u[1]+u[2]) p.L₂*cos(u[1]+u[2])   1   0;
             p.Lₐ*sin(u[1]) + p.L₂*sin(u[1]+u[2]) p.L₂*sin(u[1]+u[2])   0   1]

        M = [-m₁*(p.Lᵦ*cos(u[1]+u[2]) + (p.Lₐ-p.L₁)*cos(u[1]))  -m₁*(p.Lᵦ*sin(u[1]+u[2]) + (p.Lₐ-p.L₁)*sin(u[1]))   -m₂*((p.Lᵦ-p.L₂)*cos(u[1]+u[2]))    -m₂*((p.Lᵦ-p.L₂)*sin(u[1]+u[2]));
             -m₁*(p.Lₐ-p.L₁)*cos(u[1])                          -m₁*(p.Lₐ-p.L₁)*sin(u[1])                            0                                   0]

        C = [1 0;
             0 1;
            -p.Lᵦ*cos(u[1]+u[2]) - p.Lₐ*cos(u[1])  -p.Lᵦ*cos(u[1]+u[2]);
            -p.Lᵦ*sin(u[1]+u[2]) - p.Lₐ*sin(u[1])  -p.Lᵦ*sin(u[1]+u[2])]
        
        integrator.u[5:8] = C*inv(M*A)*M*B*u[5:8]
        integrator.p[2] = 2
    elseif idx == 2
        display("hello")
        #display(integrator.u)
        integrator.p[2] = 1
    end

end
#B = [p.L₁*cos(u[1])                       0                     1   0;
#     p.L₁*sin(u[1])                       0                     0   1;
#     p.Lₐ*cos(u[1]) + p.Lᵦ*cos(u[1]+u[2]) p.Lᵦ*cos(u[1]+u[2])   1   0;
#     p.Lₐ*sin(u[1]) + p.Lᵦ*sin(u[1]+u[2]) p.Lᵦ*sin(u[1]+u[2])   0   1]



### Convert to Cartesian ###
function X_Y(u,p)
    
     Lₐ,Lᵦ,L₁,L₂ = p
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


##################### Start of Code #####################
### Initial Conditions ###
x₀ = [-.4, -.2, .24, .5, 0, 0, 0, -10];
m₁,m₂,Lₐ,Lᵦ,L₁,L₂,c₁,c₂,g = (4,2,.2,.15,.05,.06,.2,.2,9.81)
state = 1
p = [(m₁=m₁,m₂=m₂,Lₐ=Lₐ,Lᵦ=Lᵦ,L₁=L₁,L₂=L₂,c₁=c₁,c₂=c₂, g=g) state]

####### ODE Solver #######
tspan = (0.0, 7);
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
idx = Int(dim/2+1)
plt = plot(u[1,:], u[idx,:], lc=:black, lw=2)
plot!([u[1,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₁ V.S θ̇₁")
xlabel!("θ₁")
ylabel!("θ̇₁")
display(plt)

##### θ₃ V.S θ̇₃ #####
idx = Int(dim/2+2)
plt = plot(u[2,:], u[idx,:], lc=:black, lw=2)
plot!([u[2,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot θ₃ V.S θ̇₃")
xlabel!("θ₁")
ylabel!("θ̇₃")
display(plt)

##### x V.S ẋ #####
idx = Int(dim/2+3)
plt = plot(u[3,:], u[idx,:], lc=:black, lw=2)
plot!([u[3,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot x V.S ẋ")
xlabel!("x")
ylabel!("ẋ")
display(plt)

##### y V.S ẏ #####
idx = Int(dim/2+4)
plt = plot(u[4,:], u[idx,:], lc=:black, lw=2)
plot!([u[4,1]],[u[idx,1]],markercolors = :red, shape = :circle, markersize = 3)
title!("Phase Plot y V.S ẏ")
xlabel!("y")
ylabel!("ẏ")
display(plt)


####### Create GIF #######
n =  length(sol)
p = p[1]
@gif for i ∈ 1:n
    hi = sol.u[i][4] - p.Lₐ*cos(sol.u[i][1]) - p.Lᵦ*cos(sol.u[i][1]+ sol.u[i][2])
    plot(x[:,i],y[:,i],ylim = [-.5,.5],xlim=[-.5,.5],xlabel="X Position", ylabel="Y Position", 
       markercolors = [:black, :red, :black, :red, :blue], shape = [:none,:circle, :cross, :circle, :circle],
       markersize = [3,8,6,8,16], linewidth = 4, title = round(1000*hi)/1000, aspect_ratio=:equal)
end every 1
