normsq(x) = dot(x, x)



doc"""Continued fraction algorithm: find a rational approximation $p/q$
to $x$ which is within a distance $\epsilon$ of $x$.
"""
function rational_approximation(x, ϵ)
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

doc"""
Find a collision of a particle with position `x` and velocity `v` with a local disc of radius $r$; if there is no such
collision with a nearby disc, find the nearest position of the form $(n, b)$ with $n$ an integer.
Supposes that v is in the first octant.
"""
function local_step(x, v, r)
    e1 = [1, 0]
    e2 = [0, 1]

    nn = ifloor(x)  # bottom left of current cell
    xx = x-nn       # position relative to current cell

    t = (1 - xx[1]) / v[1]  # collision time with next vertical wall
    x1 = xx + v * t         # collision position

    t2 = -xx[1] / v[1]  # collision time with previous vertical wall
    x2 = xx + v * t2

    b1 = x1[2]  # collision heights
    b2 = x2[2]

    ϵ = r / v[1]

    if (xx - e1)⋅v < 0
        if abs(b1) < ϵ
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

doc"""
Efficiently calculate with which disc of radius $r$ a trajectory will collide.
The trajectory is a line through the point $(0, b)$ with slope $m$.
"""
function efficient_disc_collision(m, b, r)
	kn = 0
    b1 = b
    ϵ = r * √(m^2 + 1)

    if b < ϵ || 1 - b < ϵ
        if b < 0.5
			q, p = rational_approximation(m, 2b)
		else
			q, p = rational_approximation(m, 2*(1 - b))
		end

		b = mod(m*q + b, 1)
		kn += q
    end

	while b > ϵ && 1 - b > ϵ
		if b < 0.5
			(q, p) = rational_approximation(m, 2b)
		else
			(q, p) = rational_approximation(m, 2*(1 - b))
		end
		b = mod(m*q + b, 1)
		kn += q
        if abs(b - b1) < 1e-30
            throw(ArgumentError("Rational velocity along channel with slope $m; no collision"))
        end
	end
	q = kn
    p = round(m*q + b1)
    return [q, p]
end

"""Find coordinates of the obstacle with which a particle
at initial position `x`, with velocity `v`, collides
if both components of the velocity are positive and the slope is less than 1,
    i.e. if the velocity vector lies within the first octant.
"""
function find_next_disc_first_octant(x, v, r)
    if normsq(round(x) - x) < r^2   # if inside obstacle
        return round(x)             # then first collision is with that obstacle
    end                             # (used for the 3D algorithm)

    x1, collided = local_step(x, v, r)
    if collided == 0  # hit local disc
        return x1
    end

    v1, v2 = v
    m = v2 / v1
    b = x1[2]
    b -= floor(b)   # move to form $(0, b)$

    centre = efficient_disc_collision(m, b, r)

    d = [0, round(b)]
    x2 = round(x1) - d + centre
    return x2
end


doc"rotation matrix clockwise by angle π/2"
const Rot = [ 0.0  1.0;
             -1.0  0.0 ]

doc"reflection matrix; maps (x,y) ↦ (y,x)"
const Refl = [ 0.0  1.0;
               1.0  0.0 ]


const T = Array(Matrix{Float64}, 8)
const T_inverse = Array(Matrix{Float64}, 8)


function initialise_transformations()
    for n in 1:8   # octant between 1 and 8

        quadrant = (n-1) ÷ 2  # integer division; quadrant between 0 and 3
        T[n] = Rot^quadrant

        if n%2 == 0
            T[n] = Refl * T[n]
        end

        T_inverse[n] = inv(T[n])
    end
end

initialise_transformations()


function octant(v)
    v1, v2 = v

    θ = atan2(v2, v1)
    if θ < 0
        θ += 2π
    end

    iceil(θ / (2π) * 8)  # from 1 to 8
end


function find_next_disc(x, v, r)

    n = octant(v)

    M = T[n]

    x = M * x
    v = M * v

    x = find_next_disc_first_octant(x, v, r)

    x = T_inverse[n] * x
end

doc"""
Find collision of particle with disc with centre `c` and radius `r`.
`x` and `v` are the position and velocity of the particle.
"""
function collision(x, c, r, v)
    v2 = normsq(v)
    b = ((x - c) ⋅ v) / v2
    c = (normsq(x - c) - r^2) / v2
    if b^2 - c < 0
        return x, +Inf   # no collision
    end

    t = -b - √(b^2 - c)  # collision time
    x = v*t + x
    return x, t

end

doc"""
Calculate the post-collision velocity at the point `x` on a disc of radius
`c`. The pre-collision velocity is `v`.
"""
function post_collision_velocity(x, c, v)
    n = x - c
    n /= norm(n)

    v -= 2(v⋅n) * n
    v /= norm(v)
    v
end

function LorentzGas(x, v, r, steps)

    for i = 1:steps
        center = find_next_disc(x, v, r)
        x, t = collision(x, center, r, v)
        v = post_collision_velocity(x, center, v)
    end

    x, v
end

function LorentzGas2(x, v, r, time)
    t = 0
    v1 = copy(v)

    while t < time
        x1 = copy(x)
        v1 = copy(v)
        center = find_next_disc(x, v, r)
        x, tt = collision(x, center, r, v)
        v = post_collision_velocity(x, center, v)
        t += tt # norm(x-x1)
    end

    t -= time
    if t==0
        v1 = v
    end
    x -= t*v1

    x, v1
end
