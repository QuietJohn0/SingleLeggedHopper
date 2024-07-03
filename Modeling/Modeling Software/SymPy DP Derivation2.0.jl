############################### Flight Phase ###############################

using SymPy

## Variables
t = Sym("t")
(θ₁, θ₂, x, y) = SymFunction("θ₁, θ₂, x, y")
(θ̇₁, θ̇₂, ẋ, ẏ) = SymFunction("θ̇₁, θ̇₂, ẋ, ẏ")
(θ̈₁, θ̈₂, ẍ, ÿ) = SymFunction("θ̈₁, θ̈₂, ẍ, ÿ")
(m₀, m₁, m₂, Lₐ, Lᵦ, L₁, L₂, I₁, I₂, c₁, c₂, cₛ₁, kₛ₁, cₛ₂, kₛ₂, g) = Sym("m₀, m₁, m₂, Lₐ, Lᵦ, L₁, L₂, I₁, I₂, c₁, c₂, cₛ₁, kₛ₁, cₛ₂, kₛ₂, g")
(τ₁, τ₂) = Sym("τ₁, τ₂")

## Time Derivitive
D(i) = diff(i,t)

## Kinetics
x₀ = x(t)
x₁ = x(t) + L₁*sin(θ₁(t))
x₂ = x(t) + Lₐ*sin(θ₁(t)) + L₂*sin(θ₁(t)+θ₂(t))

y₀ = y(t)
y₁ = y(t) - L₁*cos(θ₁(t))
y₂ = y(t) - Lₐ*cos(θ₁(t)) - L₂*cos(θ₁(t)+θ₂(t))

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
V = m₀*g*y₀ + m₁*g*y₁ + m₂*g*y₂ 

## Kinetic Energy
T = .5*m₀*(ẋ₀^2 + ẏ₀^2) + .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*m₂*(ẋ₂^2 + ẏ₂^2) + .5*I₁*θ̇cm₁^2 + .5*I₂*θ̇cm₂^2

## Lagrangian
L = T - V
w = τ₁*θ₁(t) + τ₂*θ₂(t) - c₁*θ̇₁(t)*θ₁(t) - c₂*θ̇₂(t)*θ₂(t)

## Eulur Lagrangian
Q(q) = D(diff(L,D(q))) - diff(L,q) - diff(w,q)
EL = Q.([θ₁(t); θ₂(t); x(t); y(t)])

## Subsitute
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))
EL = EL.subs(D(D(θ₂(t))), θ̈₂(t))
EL = EL.subs(D(D(x(t))), ẍ(t))
EL = EL.subs(D(D(y(t))), ÿ(t))

EL = EL.subs(D(θ₁(t)), θ̇₁(t))
EL = EL.subs(D(θ₂(t)), θ̇₂(t))
EL = EL.subs(D(x(t)), ẋ(t))
EL = EL.subs(D(y(t)), ẏ(t))

## Find State Space Equasion
M = expand.(simplify.(expand.(EL.jacobian([θ̈₁(t) θ̈₂(t) ẍ(t) ÿ(t)]))))
B = M*[θ̈₁(t); θ̈₂(t); ẍ(t); ÿ(t)]
u = [τ₁; τ₂; 0; 0]
V = simplify.(expand.(EL-B+u))

## Change the current working directory to the desired folder
parent_path = joinpath(@__DIR__, "..")
desired_folder_path = joinpath(parent_path, "Model Equasions")
cd(desired_folder_path)

## Write the equations to a text file
open("Flight.txt", "w") do file
    write(file, "")
end

global Vs = []
for i in 1:length(V)
    Vi = string(V[i,1])
    Vi = replace(Vi, "*" => "*s.", " g" => " s.g")
    Vi = replace(Vi, "1.0*" => "","2.0*" => "2*", "(t)" => "")
    Vi = replace(Vi,"θ₁" => "u[1]", "θ₂" => "u[2]", "θ̇₁" => "u[5]", "θ̇₂" => "u[6]", "x" => "u[3]", "ẋ" => "u[7]", "y" => "u[4]", "ẏ" => "u[8]")
    Vi = replace(Vi, "s.u" => "u", r"s.(\w+)\(" => s"\1(")
    
    open("Flight.txt", "a") do file
        write(file, "v" * string(i) * " = " * Vi * "\n")
    end    
    
    global Vs = push!(Vs,Vi)
end

open("Flight.txt", "a") do file
    write(file, "\n\n\n")
end  

global Mss = [;]
for i in 1:length(M[:,1])
    Mis = []
    for j in 1:length(M[1,:])
        Mij = string(M[i,j])
        Mij = replace(Mij, "*" => "*s.")
        Mij = replace(Mij, "1.0*" => "","2.0*" => "2*", "(t)" => "")
        Mij = replace(Mij,"θ₁" => "u[1]", "θ₂" => "u[2]", "θ̇₁" => "u[5]", "θ̇₂" => "u[6]", "x" => "u[3]", "ẋ" => "u[7]", "y" => "u[4]", "ẏ" => "u[8]")
        Mij = replace(Mij, "s.u" => "u", r"s.(\w+)\(" => s"\1(")
        
        open("Flight.txt", "a") do file
            write(file, "m" * string(i) * string(j) * " = " * Mij * "\n")
        end

        Mis = push!(Mis,Mij)         
    end
    global Mss = push!(Mss,Mis)
end