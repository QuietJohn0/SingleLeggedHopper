############################ Stance Phase ###############################

using SymPy

## Variables
t = Sym("t")
(θ₁, θ₂) = SymFunction("θ₁, θ₂")
(yₛ, ẏₛ, ÿₛ)  = SymFunction("yₛ, ẏₛ, ÿₛ")
(xₛ, ẋₛ, ẍₛ)  = SymFunction("xₛ, ẋₛ, ẍₛ")
(θ̇₁, θ̇₂) = SymFunction("θ̇₁, θ̇₂")
(θ̈₁, θ̈₂) = SymFunction("θ̈₁, θ̈₂")
(m₀, m₁, m₂, Lₐ, Lᵦ, L₁, L₂, I₁, I₂, c₁, c₂, cₛₓ, kₛₓ, cₛ, kₛ, g) = Sym("m₀, m₁, m₂, Lₐ, Lᵦ, L₁, L₂, I₁, I₂, c₁, c₂, cₛₓ, kₛₓ, cₛ, kₛ, g")
(τ₁, τ₂) = Sym("τ₁, τ₂")
(x_foot, y_foot) = Sym("x_foot, y_foot")

## Time Derivitive
D(i) = diff(i,t)

## Kinematics
x₀ = x_foot + xₛ(t)- Lᵦ*sin(θ₁(t)+θ₂(t)) - Lₐ*sin(θ₁(t))
x₁ = x_foot + xₛ(t) - Lᵦ*sin(θ₁(t)+θ₂(t)) - (Lₐ-L₁)*sin(θ₁(t))
x₂ = x_foot + xₛ(t)- (Lᵦ-L₂)*sin(θ₁(t)+θ₂(t))

y₀ = y_foot + yₛ(t) + Lᵦ*cos(θ₁(t)+θ₂(t)) + Lₐ*cos(θ₁(t))
y₁ = y_foot + yₛ(t) + Lᵦ*cos(θ₁(t)+θ₂(t)) + (Lₐ-L₁)*cos(θ₁(t))
y₂ = y_foot + yₛ(t) + (Lᵦ-L₂)*cos(θ₁(t)+θ₂(t))

## Velocities
ẋ₀ = D(x₀)
ẋ₁ = D(x₁)
ẋ₂ = D(x₂)
ẏ₀ = D(y₀)
ẏ₁ = D(y₁)
ẏ₂ = D(y₂)
θ̇cm₁ = D(θ₁(t))
θ̇cm₂ = D(θ₁(t))+D(θ₂(t))

## Potental Energy
V = m₀*g*y₀ + m₁*g*y₁ + m₂*g*y₂ + .5*kₛₓ*xₛ(t)^2 + .5*kₛ*yₛ(t)^2

## Kinetic Energy
T = .5*m₀*(ẋ₀^2 + ẏ₀^2) + .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*m₂*(ẋ₂^2 + ẏ₂^2)  + .5*I₁*θ̇cm₁^2 + .5*I₂*θ̇cm₂^2

## Lagrangian
L = T - V
w = τ₁*θ₁(t) + τ₂*θ₂(t) - c₁*θ̇₁(t)*θ₁(t) - c₂*θ̇₂(t)*θ₂(t) - cₛ*yₛ(t)*ẏₛ(t) - cₛₓ*xₛ(t)*ẋₛ(t)

## Eulur Lagrangian
Q(q) = D(diff(L,D(q))) - diff(L,q) - diff(w,q)
EL = Q.([θ₁(t); θ₂(t); xₛ(t); yₛ(t)])

## Subsitute
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))
EL = EL.subs(D(D(θ₂(t))), θ̈₂(t))
EL = EL.subs(D(D(xₛ(t))), ẍₛ(t))
EL = EL.subs(D(D(yₛ(t))), ÿₛ(t))

EL = EL.subs(D(θ₁(t)), θ̇₁(t))
EL = EL.subs(D(θ₂(t)), θ̇₂(t))
EL = EL.subs(D(xₛ(t)), ẋₛ(t))
EL = EL.subs(D(yₛ(t)), ẏₛ(t))

## Find State Space Equasion
M = expand.(simplify.(expand.(EL.jacobian([θ̈₁(t) θ̈₂(t) ẍₛ(t) ÿₛ(t)]))))
B = M*[θ̈₁(t); θ̈₂(t); ẍₛ(t); ÿₛ(t)]
u = [τ₁; τ₂; 0; 0]
V = simplify.(expand.(EL-B+u))

## Write new equations to file
open("eqns.txt", "w") do file
    write(file, "")
end

## Append Vs to file
global Vs = []
for i in 1:length(V)
    Vi = string(V[i,1])
    Vi = replace(Vi, "*" => "*p.", " g" => " p.g")
    Vi = replace(Vi, "1.0*" => "","2.0*" => "2*", "(t)" => "")
    Vi = replace(Vi,"θ₁" => "u[1]", "θ₂" => "u[2]", "xₛ" => "u[3]", "yₛ" => "u[4]", "θ̇₁" => "u[5]", "θ̇₂" => "u[6]", "ẋₛ" => "u[7]", "ẏₛ" => "u[8]")
    Vi = replace(Vi, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
    
    open("eqns.txt", "a") do file
        write(file, "v" * string(i) * " = " * Vi * "\n")
    end    
    global Vs = push!(Vs,Vi)
end

open("eqns.txt", "a") do file
    write(file, "\n\n\n")
end  

## Append Mss to file
global Mss = [;]
for i in 1:length(M[:,1])
    Mis = []
    for j in 1:length(M[1,:])
        Mij = string(M[i,j])
        Mij = replace(Mij, "*" => "*p.")
        Mij = replace(Mij, "1.0*" => "","2.0*" => "2*", "(t)" => "")
        Mij = replace(Mij,"θ₁" => "u[1]", "θ₂" => "u[2]", "xₛ" => "u[3]", "yₛ" => "u[4]", "θ̇₁" => "u[5]", "θ̇₂" => "u[6]", "ẋₛ" => "u[7]", "ẏₛ" => "u[8]")
        Mij = replace(Mij, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
        
        open("eqns.txt", "a") do file
            write(file, "m" * string(i) * string(j) * " = " * Mij * "\n")
        end
        Mis = push!(Mis,Mij)         
    end
    global Mss = push!(Mss,Mis)
end