# https://www.youtube.com/watch?v=uiQpwMQZBTA

f = (x-> x+1)
println(f(9)*pi)


function ode_fn!(du,u,p,t)
    th1,th2,th1d,th2d = u
    du[1] = th1
    du[2] = th2
    du[3] = th1d
    du[4] = th2d
end