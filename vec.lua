--vec.lua
--a is vector, b is scalar
local vec = {}
function vec.mult(a,b)
    return {a[1] * b, a[2] * b}
end
--a and b are vectors
function vec.add(a,b)
    return {a[1] + b[1], a[2] + b[2]}
end

function vec.magnitude(v)
    return math.sqrt(v[1] * v[1] + v[2] * v[2])
end
return vec