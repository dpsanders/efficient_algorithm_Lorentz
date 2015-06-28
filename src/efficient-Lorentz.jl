normsq(x) = dot(x, x)

"Find collision with disc with centre c"
function collision(x1, c, r, v)
    v2 = dot(v,v)
    b = dot((x1 - c), v) / v2
    c = (normsq(x1 - c) - r^2) / v2
    if b^2 - c < 0
        return false, 0
    end

    t = -b - √(b^2 - c)
    x = v*t + x1
    return x, t

end

function velo_col(x1, c, v)
    n = x1 - c
    n /= norm(n)

    v -= 2*(n⋅v)*n
    v /= norm(v)
    v
end

function frac(x, ϵ)
   h1, h2 = 1, 0
   k1, k2 = 0, 1
   b = x
   while abs(k1*x - h1) > ϵ
       a = ifloor(b)
       h1, h2 = a*h1 + h2, h1
       k1, k2 = a*k1 + k2, k1
       b = 1/(b - a)
   end
   return k1, h1
end

function nextt1(x,v,r)
    e1 = [1,0]
    e2 = [0,1]

    nn = ifloor(x)
    xx = x-nn
    t = (1 - xx[1]) / v[1]
    t2 = -xx[1]/v[1]
    x1 = xx + v*t
    x2 = xx + v*t2
    b1 = x1[2]
    b2 = x2[2]
    ϵ=r/v[1]
    testt=0

    if( dot(xx-e1,v)<0)
        if(abs(b1)<ϵ)
            return e1+nn, 0
        end
    end
    if( dot(xx-e2,v)<0)
        if(abs(1-b2)<ϵ)
            return e2+nn, 0
        end
    end

    if( dot(xx-e1-e2,v)<0)
        if(abs(1-b1)<ϵ)
            return e1+e2+nn, 0
        end
    end
    if( dot(xx-e1-e2-e2,v)<0)
        if(abs(2-b1)<ϵ)
            return e1+e2+e2+nn, 0
        end
    end
    testt=1
    return x1+nn, testt
end

function eff2(m, b, r)
	kn = 0
    b1 = copy(b)
    ϵ = r*sqrt(m*m+1)

    if(b < ϵ || 1 - b < ϵ)
        if b < 0.5
			(q, p) = frac(m, 2.*b)
		else
			(q, p) = frac(m, 2.*(1. - b))
		end
		b = mod(m*q + b, 1)
		kn += q
    end

	while b > ϵ && 1 - b > ϵ
		if b < 0.5
			(q, p) = frac(m, 2.*b)
		else
			(q, p) = frac(m, 2.*(1. - b))
		end
		b = mod(m*q + b, 1)
		kn += q
        if abs(b - b1)< 1e-30
            return false
        end
	end
	q = kn
    p = int(m*q+b1)
    return [q, p]
end

"""`Lor2` finds the coordinates of the obstacle with which a particle
at initial position `x`, with velocity `v`, collides
if both components of the velocity are positive and the slope is less than 1.
"""
function Lor2(x, v, r)
    x1, test = nextt1(x,v,r)
    if test==0
        return x1
    end

    v1 = v[1]
    v2 = v[2]
    m = v2 / v1
    b = x1[2]
    b -= floor(b)
    de=[0, int(b)]
    centre=eff2(m,b,r)

    if(centre==false)
        return false
    end

    x2=int(x1)-de+centre
    return x2
end


"rotation matrix by angle π/2"
const ROT = [ 0  1;
             -1  0 ]

"reflection matrix; interchanges y->x and x->y"

const REFL=[0  1;
           1  0 ]

const ROT2 = ROT^2
const ROT3 = ROT^3

function Lorentz2(x, v, r)
    v1, v2 = v

    vv=copy(v)
    xx=copy(x)
    m1 = v2/v1

    if norm(int(x)-x) < r
        # if a particle begins inside an obstacle, then the first collision
        # is considered with the same obstacle
        return int(x)
    end

    if m1>=0 && v2>=0  # velocity in  quadrant I
        if m1 <= 1
            xx = Lor2(xx,vv,r)
            if xx==false
                return false
            end

        elseif m1 > 1
            xx = REFL*xx
            vv = REFL*vv
            xx = Lor2(xx,vv,r)
            if xx==false
                return false
            end
            xx = REFL*xx
            vv = REFL*vv
        end
        return xx

    elseif m1 >= 0 && v2 < 0   #velocity in quadrant III
        xx = ROT2*xx
        vv = ROT2*vv
        if m1<=1
            xx = Lor2(xx,vv,r)
            if xx == false
                return false
            end

        elseif m1 > 1
            xx = REFL*xx
            vv = REFL*vv
            xx = Lor2(xx,vv,r)
            if xx == false
                return false
            end
            xx = REFL*xx
            vv = REFL*vv
        end
        xx = ROT2*xx
        vv = ROT2*vv
        return xx

    elseif m1<0 && v2>=0 # velocity in quadrant II
        xx = ROT*xx
        vv = ROT*vv
        if m1 < -1
            xx = Lor2(xx,vv,r)
            if xx == false
                return false
            end
        elseif m1 >= -1
            xx = REFL*xx
            vv = REFL*vv
            xx = Lor2(xx,vv,r)
            if xx==false
                return false
            end
            xx = REFL*xx
            vv = REFL*vv
        end
        xx = ROT3*xx
        vv = ROT3*vv
        return xx

    elseif m1<0 && v2<0  # velocity in quadrant IV
        xx = ROT3*xx
        vv = ROT3*vv
        if m1<-1
            xx = Lor2(xx,vv,r)
            if xx == false
                return false
            end

        elseif m1 >= -1
            xx = REFL*xx
            vv = REFL*vv
            xx = Lor2(xx,vv,r)
            if xx==false
                return false
            end
            xx = REFL*xx
            vv = REFL*vv
        end
        xx = ROT*xx
        vv = ROT3*vv
        return xx
    end
end

function LorentzGas1(x, v, r, steps)

    for i=1:steps
        center = Lorentz2(x, v, r)
        x, t = collision(x, center, r, v)
        v = velo_col(x, center, v)
    end

    x, v
end

function LorentzGas2(x, v, r, time)
    t=0
    v1=copy(v)

    while t<time
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

function Lorentz3D2(x, v, r)
    if norm(int(x)-x) < r   #if a particle begin inside an obstacle, then the first collision
                          #is considered with the same obstacle.
    xx, t = collision(x, int(x), r, v)
        return xx
    end

    x1=zeros(2)
    x2=zeros(2)
    x3=zeros(2)
    xr=false

    while (xr==false)
        v_xy = [v[1], v[2]]
        v_xz = [v[1], v[3]]
        v_yz = [v[2], v[3]]                #Velocities in the planes xy, xz, and yz

        vn_xy = v_xy./norm(v_xy)
        vn_xz = v_xz./norm(v_xz)
        vn_yz = v_yz./norm(v_yz,v_yz) #normalized velocities in the planes xy, xz, and yz

        tmax = 0
        tmin = -2
        x1[1] = 1
        x1[2] = 2
        x2[1] = 3
        x2[2] = 4
        x3[1] = 5
        x3[2] = 6
        #        while (int(tmax)>int(tmin+BigFloat("0.5")) )   #it can be considerin time conincidence

        while x1[1]!=x2[1] || x1[2]!=x3[1] || x2[2]!=x3[2]    #or considering the obstacle coincidence.
            x_xy=[x[1], x[2]]
            x_xz=[x[1], x[3]]
            x_yz=[x[2], x[3]]                #positions in the planes xy, xz, and yz

            x1 = Lorentz2(x_xy, vn_xy, r)#calculate the position of obstacle that collide with the particle
            if x1 == false
                return false
            end

            x2 = Lorentz2(x_xz, vn_xz, r)# in the 2D Lorentz gas for the 3 planes
            if x2 == false
                return false
            end

            x3 = Lorentz2(x_yz, vn_yz, r)
            if x3 == false
                return false
            end

            # calculate the time needed to reach the obstacle in each plane:
            t1 = (norm(x1-x_xy)-r) / norm(v_xy)
            t2 = (norm(x2-x_xz)-r) / norm(v_xz)
            t3 = (norm(x3-x_yz)-r) / norm(v_yz)

            tmax = max(t1, t2, t3)    # if the minumum and maximum valueds of time are almost the same, then
            tmin = min(t1, t2, t3)    # there is very high probable a collision
            if(tmax < 0)
                tmax = 0.1*r
            end
            x += v*(tmax-0.1*r)
        end
        xr, tt = collision(x, int(x), r, v) #calculates the collision. In case that there is not, then, it advance the particle to the
                                   #next cell
        if tt<0
            xr=false
        end
        x += v*(1-2*r)
    end
    return xr
end

function ClassicLorentz(x,v,r)    #This is a simpler and inefficient version of the Classic algorithm for ND. Just to
                                  #check if the 3D version works correctly.
   if(norm(int(x)-x)<r)   #if a particle begin inside an obstacle, then the first collision
                          #is considered with the same obstacle.
        xx,tt = collision(x,int(x),r,v)
        return xx
    end
    xr=false
    i=0.0

    while xr == false
        while xr == false
            xr, tt = collision(x,int(x+v*(i)),r,v)
            i += 1-2*r
        end

        t = norm(xr-x) / normsq(v)
        if normsq(x+t*v-xr) > normsq(xr - x)
            xr=false
        end
    end

    xr
end

function reordenar(A)
    B = zeros(0)
    N = length(A)
    for i=1:N
        ind=indmin(A)
        b=splice!(A,ind)
        push!(B,b)
    end
    return B
end

dist(x1,x2)=norm(x1-x2)
