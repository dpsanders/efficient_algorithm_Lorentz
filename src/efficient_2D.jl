normsq(x) = dot(x, x)

doc"""
Find collision of particle with disc with centre `c` and radius `r`.
`x` and `v` are the position and velocity of the particle.
"""
function collision(x, c, r, v)
    v2 = normsq(v)
    b = ((x - c) ⋅ v) / v2
    c = (normsq(x - c) - r^2) / v2
    if b^2 - c < 0  # no collision
        return false, 0
    end

    t = -b - √(b^2 - c)  # collision time
    x = v*t + x1
    return x, t

end

doc"""
Calculate the post-collision velocity at the point `x` on a disc of radius
`c`. The pre-collision velocity is `v`.
"""
function velo_col(x, c, v)
    n = x - c
    n /= norm(n)

    v -= 2(v⋅n) * n
    v /= norm(v)
    v
end

doc"""Continued fraction algorithm: find a rational approximation $p/q$
which is within a distance `ϵ` of `x`.
"""
function frac(x, ϵ)
   h1, h2 = 1, 0
   k1, k2 = 0, 1
   b = x

   while abs(k1*x - h1) > ϵ
       a = ifloor(b)
       h1, h2 = a*h1 + h2, h1
       k1, k2 = a*k1 + k2, k1
       b = inv(b - a)
   end
   return k1, h1
end


function next(x, v, r)
    e1 = [1, 0]
    e2 = [0, 1]

    nn = ifloor(x)  # bottom left of current cell
    xx = x-nn       # position relative to current cell

    t = (1 - xx[1]) / v[1]  # collision time with next vertical wall
    x1 = xx + v * t         # collision position

    t2 = -xx[1] / v[1]  # collision time with previous vertical wall
    x2 = xx + v * t2

    y1 = x1[2]  # collision heights
    y2 = x2[2]

    ϵ = r / v[1]

    if (xx - e1)⋅v < 0
        if abs(y1) < ϵ
            return e1 + nn, 0
        end
    end

    if (xx - e2)⋅v < 0
        if abs(1 - b2) <ϵ
            return e2 + nn, 0
        end
    end

    if (xx - e1 - e2)⋅v < 0
        if abs(1-b1) < ϵ
            return e1 + e2 + nn, 0
        end
    end
    if (xx - e1 - 2*e2) ⋅ v < 0
        if abs(2-b1) < ϵ
            return e1 + 2*e2 +nn, 0
        end
    end

    return x1+nn, 1
end

function eff(m, b, r)
	kn = 0
    b1 = b
    ϵ = r * √(m^2 + 1)

    if b < ϵ || 1 - b < ϵ
        if b < 0.5
			q, p = frac(m, 2b)
		else
			q, p = frac(m, 2*(1 - b))
		end

		b = mod(m*q + b, 1)
		kn += q
    end

	while b > ϵ && 1 - b > ϵ
		if b < 0.5
			(q, p) = frac(m, 2b)
		else
			(q, p) = frac(m, 2*(1 - b))
		end
		b = mod(m*q + b, 1)
		kn += q
        if abs(b - b1) < 1e-30
            throw(ArgumentError("Rational velocity along a channel: no collision"))
        end
	end
	q = kn
    p = int(m*q+b1)
    return [q, p]
end

"""Find coordinates of the obstacle with which a particle
at initial position `x`, with velocity `v`, collides
if both components of the velocity are positive and the slope is less than 1.
"""
function Lor2(x, v, r)
    x1, test = nextt1(x,v,r)
    if test==0
        return x1
    end

    v1, v2 = v
    m = v2 / v1
    b = x1[2]
    b -= floor(b)

    centre = eff2(m,b,r)

    d = [0, int(b)]
    x2 = int(x1) - d + centre
    return x2
end


"rotation matrix by angle π/2"
const ROT = [ 0.0  1.0;
             -1.0  0.0 ]

"reflection matrix; interchanges y->x and x->y"
const REFL = [ 0.0  1.0;
               1.0  0.0 ]


const T = Array(8, Array{Float,2})
const T_inverse = Array(8, Array{Float,2})


# transformations for octants:
for i in 0:7
    if i%2 == 1
        T[i] = ROT^((n-1)/2)
    else
        T[i] = REFL * ROT^(n/2)
    end

    T_inverse[i] = inv(T[i])
end


function octant(v)
    v1, v2 = v

    θ = atan2(v2, v1)
    if θ < 0
        θ += 2π
    end

    ifloor(theta / (2π) * 8)
end

function Lorentz2(x, v, r)

    n = octant(v)

    M = T[n]

    x = M * x
    v = M * v

    x = Lor2(x, v, r)

    if x == false
        return false
    end

    x = T_inverse[n] * x
end

function LorentzGas1(x, v, r, steps)

    for i = 1:steps
        center = Lorentz2(x, v, r)
        x, t = collision(x, center, r, v)
        v = velo_col(x, center, v)
    end

    x, v
end

function LorentzGas2(x, v, r, time)
    t = 0
    v1 = copy(v)

    while t < time
        x1 = copy(x)
        v1 = copy(v)
        center = Lorentz2(x, v, r)
        x, tt = collision(x, center, r, v)
        v = velo_col(x, center, v)
        t += tt # norm(x-x1)
    end

    t -= time
    if t==0
        v1 = v
    end
    x -= t*v1

    x, v1
end
