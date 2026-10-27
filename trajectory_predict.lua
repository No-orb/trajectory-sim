--trajectory_predict.lua
local vec = require('vec')
local rungekutta = require('rungekutta')
local RK4 = rungekutta.RK4
local add = vec.add
local mult = vec.mult
local magnitude = vec.magnitude
local rho = 1.225 --ground level air density in kg/m^3, ground level assumed, constant air density assumed
local trajectory_predict = {}


--p0 and v0 are a vector: [x, y], in meters
--mass in kg
--radius in meters
--h is the RK4 step size in seconds
--stop_height is the y-value at which we stop integrating (default 0, i.e. ground)
--max_steps is a safety cap so we can't loop forever (default 100000)
--ground level assumed, constant air density assumed
--returns positions, velocities, times (all parallel arrays)
function trajectory_predict.trajectory_sphere(p0, v0, mass, radius, h, stop_height, max_steps)
    assert(h > 0, "step size h must be positive")
    local c_d = .47 --Sphere drag coefficient
    local area = math.pi * radius^2
    local positions = {p0}
    local velocities = {v0}
    local times = {0}
    local stop_height = stop_height or 0
    local max_steps = max_steps or 100000

    --acceleration as a function of velocity only (autonomous system, no t dependence)
    local function acceleration(v)
        local vmag = magnitude(v)
        local dragMag = 0.5 * c_d * rho * area * vmag * vmag -- 1/2*(C_d *ρ*A*v^2), magnitude of drag force
        local drag = (vmag > 0) and mult(v, -dragMag / vmag) or {0, 0} --drag force opposes motion; avoid division by 0 at rest
        local gravity = {0, -9.807 * mass}   -- gravity force (F = mg, pointing down)
        return mult(add(drag, gravity), 1/mass)  -- F = ma => a = F/m
    end

   
    -- derivative of the pair: (p, v) -> (dp/dt, dv/dt)
    local function deriv(t, p, v)
        return v, acceleration(v)   -- dp/dt = v,  dv/dt = a(v)
    end

    local t = 0
    local i = 1
    while i <= max_steps do
        --grab the current state
        local v_i = velocities[i]
        local p_i = positions[i]

        --take one RK4 step forward
        local p_next, v_next, t_next = RK4(p_i, v_i, t, deriv, h)

        if p_next[2] <= stop_height then ---backwards interpolation to stop height if we go past it
            local p_prev = p_i
            --fraction of the step at which we crossed stop_height; guard against 0/0
            local denom = p_prev[2] - p_next[2]
            local alpha = (denom ~= 0) and (p_prev[2] - stop_height) / denom or 0
            --linear interpolation of x position at the crossing
            local x_hit = p_prev[1] + alpha * (p_next[1] - p_prev[1])
            positions[i+1]  = { x_hit, stop_height }
            --linear interpolation of velocity at the crossing
            velocities[i+1] = {
                v_i[1] + alpha * (v_next[1] - v_i[1]),
                v_i[2] + alpha * (v_next[2] - v_i[2]),
            }
            --interpolated time of impact (must be set here since we break below)
            times[i+1] = t + alpha * h
            return positions, velocities, times
        end

        --didn't hit yet: record the full step and advance
        positions[i+1]  = p_next
        velocities[i+1] = v_next
        times[i+1]      = t_next
        t = t_next
        i = i + 1
    end

    --hit the step cap without reaching stop_height: don't silently return a bad trajectory
    error("trajectory_sphere: exceeded max_steps without reaching stop_height")
end

return trajectory_predict