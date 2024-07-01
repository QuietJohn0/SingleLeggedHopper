############## Stance Phase ####################

using SymPy

## Variables
t = Sym("t")
(x, θ₁) = SymFunction("x, θ₁")
(ẋ, θ̇₁) = SymFunction("ẋ, θ̇₁")
(ẍ, θ̈₁) = SymFunction("ẍ, θ̈₁")
(M₁, m₁, I₁, La, L₁, b₁, c₁, g) = Sym("M1, m1, I1, LA, L1, b1, c1, g")
(Fₓ) = Sym("Fx")
(y₀) = Sym("yo")

s = 1
## Time Derivitive
D(i) = diff(i,t)

##Kinetics
x₁ = x(t) + L₁*sin(θ₁(t))*s
y₁ = y₀   - L₁*cos(θ₁(t))*s

## Velocities
ẋₘ = D(x(t))
ẋ₁ = D(x₁)
ẏ₁ = D(y₁)
θ̇cm₁ = D(θ₁(t))

## Potental Energy
V = m₁*g*y₁

## Kinetic Energy
T = .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*I₁*θ̇cm₁^2 + .5*M₁*ẋₘ^2

## Lagrangian
L = T - V
w = Fₓ*x(t) - b₁*x(t)*ẋ(t)- c₁*θ̇₁(t)*θ₁(t)

## Eulur Lagrangian
Q(q) = D(diff(L,D(q))) - diff(L,q) - diff(w,q)
EL = Q.([x(t); θ₁(t)])

## Subsitute
EL = EL.subs(D(D(x(t))), ẍ(t))
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))

EL = EL.subs(D(x(t)), ẋ(t))
EL = EL.subs(D(θ₁(t)), θ̇₁(t))

## Find State Space Equasion
## M*q̈ + V = F
M = expand.(simplify.(expand.(EL.jacobian([ẍ(t) θ̈₁(t)]))))
B = M*[ẍ(t); θ̈₁(t)]
F = [Fₓ; 0]
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