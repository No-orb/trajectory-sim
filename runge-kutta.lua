--runge-kutta.lua
--y is a scalar
--Returns the next y value and the next t value, dy/dt = f(t, y)
local vec = require('vec')
local add = vec.add
local mult = vec.mult
local rungekutta = {}

 function rungekutta.RK4(p, v, t, deriv, h) -- doing both the position and velocity vectors in the same step
    -- k1
    local k1p, k1v = deriv(t,         p,             v)
    -- k2
    local k2p, k2v = deriv(t + h/2,   add(p, mult(k1p, h/2)),
                                      add(v, mult(k1v, h/2)))
    -- k3
    local k3p, k3v = deriv(t + h/2,   add(p, mult(k2p, h/2)),
                                      add(v, mult(k2v, h/2)))
    -- k4
    local k4p, k4v = deriv(t + h,     add(p, mult(k3p, h)),
                                      add(v, mult(k3v, h)))

    local p_next = add(p, mult(add(k1p, add(mult(k2p, 2), add(mult(k3p, 2), k4p))), h/6))
    local v_next = add(v, mult(add(k1v, add(mult(k2v, 2), add(mult(k3v, 2), k4v))), h/6))
    return p_next, v_next, t + h
end

return rungekutta