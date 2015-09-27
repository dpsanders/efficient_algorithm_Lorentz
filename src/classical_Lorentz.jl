function ClassicLorentz(x, v, r)    #This is a simpler and inefficient version of the Classic algorithm for ND. Just to
                                  #check if the 3D version works correctly.
   if norm(int(x) - x) < r    #if a particle begin inside an obstacle, then the first collision
                          #is considered with the same obstacle.
        xx,tt = collision(x, int(x), r, v)
        return xx
    end
    xr = false
    i = 0.0

    while xr == false
        while xr == false
            xr, tt = collision(x, int(x+v*(i)), r, v)
            i += 1 - 2*r
        end

        t = norm(xr - x) / normsq(v)
        if normsq(x + t*v - xr) > normsq(xr - x)
            xr = false
        end
    end

    xr
end
