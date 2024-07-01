using SymPy

## Variables
t = Sym("t")
(θ₁, θ₃, x, y) = SymFunction("θ₁, θ₃, x, y")
(θ̇₁, θ̇₃, ẋ, ẏ) = SymFunction("θ̇₁, θ̇₃, ẋ, ẏ")
(θ̈₁, θ̈₃, ẍ, ÿ) = SymFunction("θ̈₁, θ̈₃, ẍ, ÿ")
(m₁, m₂, Lₐ, Lᵦ, L₁, L₂, c₁, c₂, kₛ, lₛ, g) = Sym("m₁, m₂, Lₐ, Lᵦ, L₁, L₂, c₁, c₂, kₛ, lₛ, g")
(τ₁, τ₃) = Sym("τ₁, τ₃")

## Time Derivitive
D(i) = diff(i,t)

## Kinetics
x₁ = x(t) + L₁*sin(θ₁(t))
x₂ = x(t) + Lₐ*sin(θ₁(t)) + L₂*sin(θ₁(t)+θ₃(t))

y₁ = y(t) - L₁*cos(θ₁(t))
y₂ = y(t) - Lₐ*cos(θ₁(t)) - L₂*cos(θ₁(t)+θ₃(t))

## Velocities
ẋ₁ = D(x₁)
ẋ₂ = D(x₂)
ẏ₁ = D(y₁)
ẏ₂ = D(y₂)

## Potental Energy
V = m₁*g*y₁ + m₂*g*y₂

## Kinetic Energy
T = .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*m₂*(ẋ₂^2 + ẏ₂^2)

## Lagrangian
L = T - V
w = τ₁*θ₁(t) + τ₃*θ₃(t) - c₁*θ̇₁(t)*θ₁(t) - c₂*θ̇₃(t)*θ₃(t)

## Eulur Lagrangian
Q(q) = D(diff(L,D(q))) - diff(L,q) - diff(w,q)
EL = Q.([θ₁(t); θ₃(t); x(t); y(t)])

## Subsitute
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))
EL = EL.subs(D(D(θ₃(t))), θ̈₃(t))
EL = EL.subs(D(D(x(t))), ẍ(t))
EL = EL.subs(D(D(y(t))), ÿ(t))

EL = EL.subs(D(θ₁(t)), θ̇₁(t))
EL = EL.subs(D(θ₃(t)), θ̇₃(t))
EL = EL.subs(D(x(t)), ẋ(t))
EL = EL.subs(D(y(t)), ẏ(t))

## Find State Space Equasion
M = expand.(simplify.(expand.(EL.jacobian([θ̈₁(t) θ̈₃(t) ẍ(t) ÿ(t)]))))
B = M*[θ̈₁(t); θ̈₃(t); ẍ(t); ÿ(t)]
τ = [τ₁; τ₃; 0; 0]
V = simplify.(expand.(EL-B+τ))

open("eqns.txt", "w") do file
    write(file, "")
end

global Vs = []
for i in 1:length(V)
    Vi = string(V[i,1])
    Vi = replace(Vi, "*" => "*s.", " g" => " s.g")
    Vi = replace(Vi, "1.0*" => "","2.0*" => "2*", "(t)" => "")
    Vi = replace(Vi,"θ₁" => "u[1]", "θ₃" => "u[2]", "θ̇₁" => "u[5]", "θ̇₃" => "u[6]", "x" => "u[3]", "ẋ" => "u[7]", "y" => "u[4]", "ẏ" => "u[8]")
    Vi = replace(Vi, "s.u" => "u", r"s.(\w+)\(" => s"\1(")
    
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
        Mij = replace(Mij, "*" => "*s.")
        Mij = replace(Mij, "1.0*" => "","2.0*" => "2*", "(t)" => "")
        Mij = replace(Mij,"θ₁" => "u[1]", "θ₃" => "u[2]", "θ̇₁" => "u[5]", "θ̇₃" => "u[6]", "x" => "u[3]", "ẋ" => "u[7]", "y" => "u[4]", "ẏ" => "u[8]")
        Mij = replace(Mij, "s.u" => "u", r"s.(\w+)\(" => s"\1(")
        
        open("eqns.txt", "a") do file
            write(file, "m" * string(i) * string(j) * " = " * Mij * "\n")
        end

        Mis = push!(Mis,Mij)         
    end
    global Mss = push!(Mss,Mis)
end