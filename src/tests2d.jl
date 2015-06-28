#include("efficient-algorithm-revised.jl")
include("eff-ata-latest-old.jl")
push!(LOAD_PATH, "modules")

using ClassicalLorentz
using EfficientLorentz

# Ata's E vs my E
#=
for j = 10:50
    for i = 1:100
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand()
	    v = [cos(phi), sin(phi)]
        r = 0.8^j

        effic = Lorentz(x, v, r)
        class = int(ClassicLorentz(x, v, r))

        my_effic = collisions(x[1], x[2], v[1], v[2], r, 1, 64)[1][1]

        if effic == my_effic
            #println("Passed with x = $x; v = $v; r = $r")
        else
            t = norm(my_effic - x)/norm(v)
            my_class = collisions_classical(x, v, r, t + 0.1, 64)[2]
            println("Failed with x = $x; v = $v; r = $r : E = $effic, MyE = $my_effic, C = $class, MyC = $my_class")
        end
    end
end
=#

# Ata's E vs Ata's C

function compare()

    for j = 10:30 #50
        @show j
        for i = 1:100
            x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
            x = big(x)

            phi = big(2pi*rand())
            v = [cos(phi), sin(phi)]

            r = big(0.8)^j
            x

            effic = Lorentz2(x, v, r)
            class = int(ClassicLorentz(x, v, r))

            if effic == class
                #println("Passed with x = $x; v = $v; r = $r")
            else
                my_effic = collisions(x[1], x[2], v[1], v[2], r, 1, 64)[1][1]
                println("Failed with x = $x; v = $v; r = $r : E = $effic, C = $class, myE = $my_effic")
            end
        end
    end
end

@time compare()

