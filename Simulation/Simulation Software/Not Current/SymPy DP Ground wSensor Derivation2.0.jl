############## Stance Phase ####################

using SymPy

## Variables
t = Sym("t")
(θ₁, θ₃) = SymFunction("θ₁, θ₃")
(yₛ, ẏₛ, ÿₛ)  = SymFunction("yₛ, ẏₛ, ÿₛ")
(θ̇₁, θ̇₃) = SymFunction("θ̇₁, θ̇₃")
(θ̈₁, θ̈₃) = SymFunction("θ̈₁, θ̈₃")
(m₀, m₁, m₂, Lₐ, Lᵦ, L₁, L₂, I₁, I₂, c₁, c₂, cₛ, kₛ, lₛ, g) = Sym("m₀, m₁, m₂, Lₐ, Lᵦ, L₁, L₂, I₁, I₂, c₁, c₂, cₛ, kₛ, lₛ, g")
(τ₁, τ₃) = Sym("τ₁, τ₃")
(x₀, y₀) = Sym("x₀, y₀")

## Time Derivitive
D(i) = diff(i,t)

##Kinetics
x₁ = x₀ - Lᵦ*sin(θ₁(t)+θ₃(t)) - Lₐ*sin(θ₁(t))
x₁ = x₀ - Lᵦ*sin(θ₁(t)+θ₃(t)) - (Lₐ-L₁)*sin(θ₁(t))
x₂ = x₀ - (Lᵦ-L₂)*sin(θ₁(t)+θ₃(t))

y₁ = y₀ + lₛ + yₛ(t) + Lᵦ*cos(θ₁(t)+θ₃(t)) + Lₐ*cos(θ₁(t))
y₁ = y₀ + lₛ + yₛ(t) + Lᵦ*cos(θ₁(t)+θ₃(t)) + (Lₐ-L₁)*cos(θ₁(t))
y₂ = y₀ + lₛ + yₛ(t) + (Lᵦ-L₂)*cos(θ₁(t)+θ₃(t))

## Velocities
ẋ₀ = D(x₀)
ẋ₁ = D(x₁)
ẋ₂ = D(x₂)
ẏ₀ = D(y₀)
ẏ₁ = D(y₁)
ẏ₂ = D(y₂)
θ̇cm₁ = D(θ₁(t))
θ̇cm₂ = D(θ₁(t))+D(θ₃(t))

## Potental Energy
V = m₀*g*y₀ + m₁*g*y₁ + m₂*g*y₂ + .5*kₛ*yₛ(t)^2

## Kinetic Energy
T = .5*m₀*(ẋ₀^2 + ẏ₀^2) + .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*m₂*(ẋ₂^2 + ẏ₂^2)  + .5*I₁*θ̇cm₁^2 + .5*I₂*θ̇cm₂^2

## Lagrangian
L = T - V
w = τ₁*θ₁(t) + τ₃*θ₃(t) - c₁*θ̇₁(t)*θ₁(t) - c₂*θ̇₃(t)*θ₃(t) - cₛ*yₛ(t)*ẏₛ(t)

## Eulur Lagrangian
Q(q) = D(diff(L,D(q))) - diff(L,q) - diff(w,q)
EL = Q.([θ₁(t); θ₃(t); yₛ(t)])

## Subsitute
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))
EL = EL.subs(D(D(θ₃(t))), θ̈₃(t))
EL = EL.subs(D(D(yₛ(t))), ÿₛ(t))

EL = EL.subs(D(θ₁(t)), θ̇₁(t))
EL = EL.subs(D(θ₃(t)), θ̇₃(t))
EL = EL.subs(D(yₛ(t)), ẏₛ(t))

## Find State Space Equasion
M = expand.(simplify.(expand.(EL.jacobian([θ̈₁(t) θ̈₃(t) ÿₛ(t)]))))
B = M*[θ̈₁(t); θ̈₃(t); ÿₛ(t)]
τ = [τ₁; τ₃; 0]
V = simplify.(expand.(EL-B+τ))


open("eqns.txt", "w") do file
    write(file, "")
end

global Vs = []
for i in 1:length(V)
    Vi = string(V[i,1])
    Vi = replace(Vi, "*" => "*p.", " g" => " p.g")
    Vi = replace(Vi, "1.0*" => "","2.0*" => "2*", "(t)" => "")
    Vi = replace(Vi,"θ₁" => "u[1]", "θ₃" => "u[2]", "θ̇₁" => "u[4]", "θ̇₃" => "u[5]", "yₛ" => "u[3]", "ẏₛ" => "u[6]")
    Vi = replace(Vi, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
    
    open("eqns.txt", "a") do file
        write(file, "v" * string(i) * " = " * Vi * "\n")
    end    
    
    global Vs = push!(Vs,Vi)
end

open("eqns.txt", "a") do file
    write(file, "\n\n\n")
end  

global Mss = [;]
for i in 1:length(M[:,1])
    Mis = []
    for j in 1:length(M[1,:])
        Mij = string(M[i,j])
        Mij = replace(Mij, "*" => "*p.")
        Mij = replace(Mij, "1.0*" => "","2.0*" => "2*", "(t)" => "")
        Mij = replace(Mij,"θ₁" => "u[1]", "θ₃" => "u[2]", "θ̇₁" => "u[4]", "θ̇₃" => "u[5]", "yₛ" => "u[3]", "ẏₛ" => "u[6]")
        Mij = replace(Mij, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
        
        open("eqns.txt", "a") do file
            write(file, "m" * string(i) * string(j) * " = " * Mij * "\n")
        end

        Mis = push!(Mis,Mij)         
    end
    global Mss = push!(Mss,Mis)
end



#=
T = T.subs(D(θ₁(t)), θ̇₁(t))
T = T.subs(D(θ₃(t)), θ̇₃(t))
T = simplify.(expand.(T))
T = collect(T,[θ̇₁(t)*θ̇₃(t),θ̇₃(t)^2,θ̇₁(t)^2])

EL = simplify.(expand.(EL))
EL1 = collect.(EL[1],[θ̈₁(t),θ̈₃(t),Lᵦ*L₂*m₁*θ̇₃(t)*sin(θ₃(t)), Lᵦ*Lₐ*m₁*θ̇₃(t)*sin(θ₃(t))])
EL2 = collect.(EL[2],[θ̈₁(t),θ̈₃(t),Lᵦ*L₂*m₁*θ̇₃(t)*sin(θ₃(t)), Lᵦ*Lₐ*m₁*θ̇₃(t)*sin(θ₃(t))])
V = collect.(V,[Lᵦ*L₂*m₁*θ̇₃(t)*sin(θ₃(t)), Lᵦ*Lₐ*m₁*θ̇₃(t)*sin(θ₃(t)),0])

(u1, u2, u5, u6) = Sym("u1, u2, u5, u6")
V = V.subs.(θ₁(t),u1)
V = V.subs.(θ₃(t),u2)
V = V.subs.(θ̇₁(t),u5)
V = V.subs.(θ̇₃(t),u6)

M = M.subs.(θ₁(t),u1)
M = M.subs.(θ₃(t),u2)
M = M.subs.(θ̇₁(t),u5)
M = M.subs.(θ̇₃(t),u6)=#
