function Lorentz3D2(x, v, r)
    if norm(int(x) - x) < r     # if a particle begin inside an obstacle, then the first collision
                                # is considered with the same obstacle.
    xx, t = collision(x, int(x), r, v)
        return xx
    end

    x1 = zeros(2)
    x2 = zeros(2)
    x3 = zeros(2)
    xr = false

    while (xr == false)
        v_xy = [v[1], v[2]]
        v_xz = [v[1], v[3]]
        v_yz = [v[2], v[3]]                #Velocities in the planes xy, xz, and yz

        vn_xy = v_xy ./ norm(v_xy)
        vn_xz = v_xz ./ norm(v_xz)
        vn_yz = v_yz ./ norm(v_yz) #normalized velocities in the planes xy, xz, and yz

        tmax = 0
        tmin = -2
        x1[1] = 1
        x1[2] = 2
        x2[1] = 3
        x2[2] = 4
        x3[1] = 5
        x3[2] = 6
        #        while (int(tmax)>int(tmin+BigFloat("0.5")) )   #it can be considerin time conincidence

        while x1[1] != x2[1] || x1[2] != x3[1] || x2[2] != x3[2]    #or considering the obstacle coincidence.
            x_xy = [x[1], x[2]]
            x_xz = [x[1], x[3]]
            x_yz = [x[2], x[3]]                #positions in the planes xy, xz, and yz

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
            t1 = (norm(x1 - x_xy) - r) / norm(v_xy)
            t2 = (norm(x2 - x_xz) - r) / norm(v_xz)
            t3 = (norm(x3 - x_yz) - r) / norm(v_yz)

            tmax = max(t1, t2, t3)    # if the minumum and maximum valueds of time are almost the same, then
            tmin = min(t1, t2, t3)    # there is very high probable a collision
            if tmax < 0
                tmax = 0.1 * r
            end
            x += v * (tmax - 0.1*r)
        end

        xr, tt = collision(x, int(x), r, v) #calculates the collision. In case that there is not, then, it advance the particle to the
                                   #next cell
        if tt<0
            xr=false
        end
        x += v * (1 - 2*r)
    end
    return xr
end
