using SymPy

## Variables
t = Sym("t")
(θ₁, θ₃, x, y) = SymFunction("θ₁, θ₃, x, y")
(θ̇₁, θ̇₃, ẋ, ẏ) = SymFunction("θ̇₁, θ̇₃, ẋ, ẏ")
(θ̈₁, θ̈₃, ẍ, ÿ) = SymFunction("θ̈₁, θ̈₃, ẍ, ÿ")
(m₁, m₂, Lₐ, Lᵦ, L₁, L₂, c₁, c₂, cₓ, c🥰, g) = Sym("m₁, m₂, Lₐ, Lᵦ, L₁, L₂, c₁, c₂, cₓ, c🥰, g")
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
A = simplify.(expand.(EL.jacobian([θ̈₁(t) θ̈₃(t) ẍ(t) ÿ(t)])))
B = A*[θ̈₁(t); θ̈₃(t); ẍ(t); ÿ(t)]
F = simplify.(expand.(B-EL))