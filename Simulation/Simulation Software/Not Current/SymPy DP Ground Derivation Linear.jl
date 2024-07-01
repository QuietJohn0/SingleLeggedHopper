############## Stance Phase ####################

using SymPy

## Variables
t = Sym("t")
(θ₁, θ₃) = SymFunction("θ₁, θ₃")
(θ̇₁, θ̇₃) = SymFunction("θ̇₁, θ̇₃")
(θ̈₁, θ̈₃) = SymFunction("θ̈₁, θ̈₃")
(m₁, m₂, Lₐ, Lᵦ, L₁, L₂, c₁, c₂, g) = Sym("m₁, m₂, Lₐ, Lᵦ, L₁, L₂, c₁, c₂, g")
(τ₁, τ₃) = Sym("τ₁, τ₃")
(x₀, y₀) = Sym("x₀, y₀")

## Time Derivitive
D(i) = diff(i,t)

##Kinetics
x₁ = x₀ - Lᵦ*sin(θ₁(t)+θ₃(t)) - (Lₐ-L₂)*sin(θ₁(t))
x₂ = x₀ - (Lᵦ-L₂)*sin(θ₁(t)+θ₃(t))

y₁ = y₀ + Lᵦ*cos(θ₁(t)+θ₃(t)) + (Lₐ-L₂)*cos(θ₁(t))
y₂ = y₀ + (Lᵦ-L₂)*cos(θ₁(t)+θ₃(t))

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
EL = Q.([θ₁(t); θ₃(t)])

## Subsitute
EL = EL.subs(D(D(θ₁(t))), θ̈₁(t))
EL = EL.subs(D(D(θ₃(t))), θ̈₃(t))

EL = EL.subs(D(θ₁(t)), θ̇₁(t))
EL = EL.subs(D(θ₃(t)), θ̇₃(t))

## Find State Space Equasion
A = simplify.(expand.(EL.jacobian([θ̈₁(t) θ̈₃(t)])))
B = A*[θ̈₁(t); θ̈₃(t)]
F = simplify.(expand.(B-EL))











T = T.subs(D(θ₁(t)), θ̇₁(t))
T = T.subs(D(θ₃(t)), θ̇₃(t))
T = simplify.(expand.(T))
T = collect(T,[θ̇₁(t)*θ̇₃(t),θ̇₃(t)^2,θ̇₁(t)^2])

EL = simplify.(expand.(EL))
EL1 = collect.(EL[1],[θ̈₁(t),θ̈₃(t),Lᵦ*L₂*m₁*θ̇₃(t)*sin(θ₃(t)), Lᵦ*Lₐ*m₁*θ̇₃(t)*sin(θ₃(t))])
EL2 = collect.(EL[2],[θ̈₁(t),θ̈₃(t),Lᵦ*L₂*m₁*θ̇₃(t)*sin(θ₃(t)), Lᵦ*Lₐ*m₁*θ̇₃(t)*sin(θ₃(t))])

F = collect.(F,[Lᵦ*L₂*m₁*θ̇₃(t)*sin(θ₃(t)), Lᵦ*Lₐ*m₁*θ̇₃(t)*sin(θ₃(t))])

(u1, u2, u3, u4) = Sym("u1, u2, u3, u4")
F = F.subs.(sin(θ₁(t)),u1)
F = F.subs.(sin(θ₃(t)),u2)
F = F.subs.(cos(θ₁(t)),1)
F = F.subs.(cos(θ₃(t)),1)
F = F.subs.(θ̇₁(t),u3)
F = F.subs.(θ̇₃(t),u4)
F = F.subs.(sin(θ₁(t)+θ₃(t)),(u1 + u2))
F = F.subs.(cos(θ₁(t)+θ₃(t)),1)

A = A.subs.(sin(θ₁(t)),u1)
A = A.subs.(sin(θ₃(t)),u2)
A = A.subs.(cos(θ₁(t)),1)
A = A.subs.(cos(θ₃(t)),1)
A = A.subs.(θ̇₁(t),u3)
A = A.subs.(θ̇₃(t),u4)
A = A.subs.(sin(θ₁(t)+θ₃(t)),(u1 + u2))
A = A.subs.(cos(θ₁(t)+θ₃(t)),1)