############## Stance Phase ####################

using SymPy

## Variables
t = Sym("t")
(x, θ₁, θ₂, θ₃) = SymFunction("x, θ₁, θ₂, θ₃")
(ẋ, θ̇₁, θ̇₂, θ̇₃) = SymFunction("ẋ, θ̇₁, θ̇₂, θ̇₃")
(ẍ, θ̈₁, θ̈₂, θ̈₃) = SymFunction("ẍ, θ̈₁, θ̈₂, θ̈₃")
(M₁, m₁, m₂, m₃, I₁, I₂, I₃, La, Lb, Lc, L₁, L₂, L₃, b₁, c₁, c₂, c₃, g) = Sym("M1, m1, m2, m3, I1, I2, I3, LA, LB, LC, L1, L2, L3, b1, c1, c2, c3, g")
(Fₓ) = Sym("Fx")
(y₀) = Sym("yo")

## Time Derivitive
D(i) = diff(i,t)

##Kinetics
x₁ = x(t) + L₁*sin(θ₁(t))
x₂ = x(t) + La*sin(θ₁(t)) + L₂*sin(θ₁(t)+θ₂(t))
x₃ = x(t) + La*sin(θ₁(t)) + Lb*sin(θ₁(t)+θ₂(t)) + L₃*sin(θ₁(t)+θ₂(t)+θ₃(t))

y₁ = y₀   - L₁*cos(θ₁(t))
y₂ = y₀   - La*cos(θ₁(t)) - L₂*cos(θ₁(t)+θ₂(t))
y₃ = y₀   - La*cos(θ₁(t)) - Lb*cos(θ₁(t)+θ₂(t)) - L₃*cos(θ₁(t)+θ₂(t)+θ₃(t))

## Velocities
ẋₘ = D(x(t))
ẋ₁ = D(x₁)
ẋ₂ = D(x₂)
ẋ₃ = D(x₃)
ẏ₁ = D(y₁)
ẏ₂ = D(y₂)
ẏ₃ = D(y₃)
θ̇cm₁ = D(θ₁(t))
θ̇cm₂ = D(θ₁(t))+D(θ₂(t))
θ̇cm₃ = D(θ₁(t))+D(θ₂(t))+D(θ₃(t))


## Potental Energy
V = m₁*g*y₁ + m₂*g*y₂ + m₃*g*y₃

## Kinetic Energy
T =  .5*M₁*ẋₘ^2 + .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*I₁*θ̇cm₁^2 + .5*m₂*(ẋ₂^2 + ẏ₂^2) + .5*I₂*θ̇cm₂^2 + .5*m₃*(ẋ₃^2 + ẏ₃^2) + .5*I₃*θ̇cm₃^2

## Lagrangian
L = T - V
w = Fₓ*x(t) - b₁*x(t)*ẋ(t) - c₁*θ̇₁(t)*θ₁(t) - c₂*θ̇₂(t)*θ₂(t) - c₃*θ̇₃(t)*θ₃(t)

## Eulur Lagrangian
Q(q) = D(diff(L,D(q))) - diff(L,q) - diff(w,q)
EL = Q.([x(t); θ₁(t); θ₂(t); θ₃(t)])

## Subsitute
EL = EL.subs(D(D(x(t))), ẍ(t))
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))
EL = EL.subs(D(D(θ₂(t))), θ̈₂(t))
EL = EL.subs(D(D(θ₃(t))), θ̈₃(t))

EL = EL.subs(D(x(t)), ẋ(t))
EL = EL.subs(D(θ₁(t)), θ̇₁(t))
EL = EL.subs(D(θ₂(t)), θ̇₂(t))
EL = EL.subs(D(θ₃(t)), θ̇₃(t))

## Find State Space Equasion
## M*q̈ + V = F
M = expand.(simplify.(expand.(EL.jacobian([ẍ(t) θ̈₁(t) θ̈₂(t) θ̈₃(t)]))))
B = M*[ẍ(t); θ̈₁(t); θ̈₂(t); θ̈₃(t)]
F = [Fₓ; 0; 0; 0]
V = simplify.(expand.(EL-B+F))


open("eqns.txt", "w") do file
    write(file, "")
end

global Vs = []
for i in 1:length(V)
    Vi = string(V[i,1])
#    Vi = replace(Vi, "*" => "*p.", " g" => " p.g")
    Vi = replace(Vi, "1.0*" => "","2.0*" => "2*", "(t)" => "")
#    Vi = replace(Vi,"p.x" => "x", "p.θ₁" => "th1", "p.θ₂" => "th2", "p.θ₃" => "th3", "p.ẋ" => "xd", "p.θ̇₁" => "thd1", "p.θ̇₂" => "thd2", "p.θ̇₃" => "thd3")
#    Vi = replace(Vi,"x" => "x", "θ₁" => "th1", "θ₂" => "th2", "θ₃" => "th3")
#    Vi = replace(Vi, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
    
    open("eqns.txt", "a") do file
        write(file, "v" * string(i) * " = " * Vi * ";\n")
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
#        Mij = replace(Mij, "*" => "*p.")
        Mij = replace(Mij, "1.0*" => "","2.0*" => "2*", "(t)" => "")
#        Mij = replace(Mij,"p.x" => "x", "p.θ₁" => "th1", "p.θ₂" => "th2", "p.θ₃" => "th3", "p.ẋ" => "xd", "p.θ̇₁" => "thd1", "p.θ̇₂" => "thd2", "p.θ̇₃" => "thd3")
#        Mij = replace(Mij,"x" => "x", "θ₁" => "th1", "θ₂" => "th2", "θ₃" => "th3")
#        Mij = replace(Mij, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
        
        open("eqns.txt", "a") do file
            write(file, "m" * string(i) * string(j) * " = " * Mij * ";\n")
        end

        Mis = push!(Mis,Mij)         
    end
    global Mss = push!(Mss,Mis)
end








open("eqns1.txt", "w") do file
    write(file, "")
end

global Vs = []
for i in 1:length(V)
    Vi = string(V[i,1])
    Vi = replace(Vi, "*" => "*p.", " g" => " p.g")
    Vi = replace(Vi, "1.0*" => "","2.0*" => "2*", "(t)" => "")
    Vi = replace(Vi,"p.x" => "x", "p.θ₁" => "th1", "p.θ₂" => "th2", "p.θ₃" => "th3", "p.ẋ" => "xd", "p.θ̇₁" => "thd1", "p.θ̇₂" => "thd2", "p.θ̇₃" => "thd3")
    Vi = replace(Vi,"x" => "x", "θ₁" => "th1", "θ₂" => "th2", "θ₃" => "th3")
    Vi = replace(Vi, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
    
    open("eqns1.txt", "a") do file
        write(file, "v" * string(i) * " = " * Vi * ";\n")
    end    
    
    global Vs = push!(Vs,Vi)
end

open("eqns1.txt", "a") do file
    write(file, "\n\n\n")
end  

global Mss = [;]
for i in 1:length(M[:,1])
    Mis = []
    for j in 1:length(M[1,:])
        Mij = string(M[i,j])
        Mij = replace(Mij, "*" => "*p.")
        Mij = replace(Mij, "1.0*" => "","2.0*" => "2*", "(t)" => "")
        Mij = replace(Mij,"p.x" => "x", "p.θ₁" => "th1", "p.θ₂" => "th2", "p.θ₃" => "th3", "p.ẋ" => "xd", "p.θ̇₁" => "thd1", "p.θ̇₂" => "thd2", "p.θ̇₃" => "thd3")
        Mij = replace(Mij,"x" => "x", "θ₁" => "th1", "θ₂" => "th2", "θ₃" => "th3")
        Mij = replace(Mij, "p.u" => "u", r"p.(\w+)\(" => s"\1(")
        
        open("eqns1.txt", "a") do file
            write(file, "m" * string(i) * string(j) * " = " * Mij * ";\n")
        end

        Mis = push!(Mis,Mij)
    end
    global Mss = push!(Mss,Mis)
end