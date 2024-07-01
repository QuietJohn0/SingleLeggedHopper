using Symbolics

@variables t
D = Differential(t)
@variables θ₁(t) θ₃(t) x(t) y(t) 
@variables θ̇₁(t) θ̇₃(t) ẋ(t) ẏ(t) 
@variables θ̈₁(t) θ̈₃(t) ẍ(t) ÿ(t)
@variables m₁ m₂ Lₐ Lᵦ L₁ L₂ c₁ c₂ cₓ c🥰 g
@variables τ₁ τ₃

##Kinetics
x₁ = x + L₁*sin(θ₁)
x₂ = x + Lₐ*sin(θ₁) + L₂*sin(θ₁+θ₃)

y₁ = y - L₁*cos(θ₁)
y₂ = y - Lₐ*cos(θ₁) - L₂*cos(θ₁+θ₃)

ẋ₁ = expand_derivatives(D(x₁))
ẋ₂ = expand_derivatives(D(x₂))
ẏ₁ = expand_derivatives(D(y₁))
ẏ₂ = expand_derivatives(D(y₂))

## Potental Energy
V = m₁*g*y₁ + m₂*g*y₂

## Kinetic Energy
T = .5*m₁*(ẋ₁^2 + ẏ₁^2) + .5*m₂*(ẋ₂^2 + ẏ₂^2)
#T = simplify(substitute(T,(Dict(D(θ₁) => θ̇₁, D(θ₃) => θ̇₃, D(x) => ẋ, D(y) => ẏ))))

## Lagrangian
L = T - V
w = τ₁*θ₁ + τ₃*θ₃ - c₁*θ̇₁*θ₁ - c₂*θ̇₃*θ₃

## Eulur Lagrangian
EL = convert(Array{Num},zeros(4,1))
## q = θ₁
q = θ₁
Q = Differential(q)
Q̇ = Differential(D(q))
Lq = D(Q̇(L)) - Q(L) - Q(w)
EL[1] = substitute(expand_derivatives(expand(Lq)),(Dict(D(D(θ₁)) => θ̈₁, D(D(θ₃)) => θ̈₃, D(D(x)) => ẍ, D(D(y)) => ÿ, D(θ₁) => θ̇₁, D(θ₃) => θ̇₃, D(x) => ẋ, D(y) => ẏ)))

## q = θ₃
q = θ₃
Q = Differential(q)
Q̇ = Differential(D(q))
Lq = D(Q̇(L)) - Q(L) - Q(w)
EL[2] = substitute(expand_derivatives(expand(Lq)),(Dict(D(D(θ₁)) => θ̈₁, D(D(θ₃)) => θ̈₃, D(D(x)) => ẍ, D(D(y)) => ÿ, D(θ₁) => θ̇₁, D(θ₃) => θ̇₃, D(x) => ẋ, D(y) => ẏ)))


## q = x
q = x
Q = Differential(q)
Q̇ = Differential(D(q))
Lq = D(Q̇(L)) - Q(L) - Q(w)
EL[3] = substitute(expand_derivatives(expand(Lq)),(Dict(D(D(θ₁)) => θ̈₁, D(D(θ₃)) => θ̈₃, D(D(x)) => ẍ, D(D(y)) => ÿ, D(θ₁) => θ̇₁, D(θ₃) => θ̇₃, D(x) => ẋ, D(y) => ẏ)))

## q = y
q = y
Q = Differential(q)
Q̇ = Differential(D(q))
Lq = D(Q̇(L)) - Q(L) - Q(w)
EL[4] = substitute(expand_derivatives(expand(Lq)),(Dict(D(D(θ₁)) => θ̈₁, D(D(θ₃)) => θ̈₃, D(D(x)) => ẍ, D(D(y)) => ÿ, D(θ₁) => θ̇₁, D(θ₃) => θ̇₃, D(x) => ẋ, D(y) => ẏ)))



A = Symbolics.expand(Symbolics.jacobian(EL, [θ̈₁ θ̈₃ ẍ ÿ]))

B = A*[θ̈₁ θ̈₃ ẍ ÿ]'
B = Symbolics.expand(B)

F = simplify(collect(expand(B - EL)))