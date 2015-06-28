module ClassicalLorentz
export crossz, crossing, crossing3d, collisions_classical, first_collision_classical, collisions3d_classical, first_collision3d_classical

# Perpendicular (z-) component of cross product (scalar quantity) for 2-D vectors: (x cross y) dot e_z
#set_bigfloat_precision(1024)
crossz(x, y) = x[1]*y[2] - x[2]*y[1]

function crossing(x, v, n, m)
	# Outputs new position in square [-0.5, 0.5)^2 and updates the numbers of new square
	
	# Minimum positive times to nearest crossing - concise formulas p. 48 - 49
	array_times = BigFloat[]
	k = v[2]/v[1]	
	for i = 1:2
		for j = 1:2
			push!(array_times, ((-1)^j*0.5 - x[i])/v[i])
		end	
	end
	
	# Extract the minimum positive time value
	minpos = Inf
	number = 0
	for i = 1:length(array_times)
		if array_times[i] < minpos && array_times[i] > 0
			minpos = array_times[i]
			number = i
		end				
	end
	
	t = array_times[number]
	x1 = x + v*t
	d = norm(v*t)
	
	# Possible outcomes depending on which wall will be hit
	if number == 2
		n += 1
		return [-0.5, x1[2]], d, n, m
	elseif number == 4
		m += 1
		return [x1[1], -0.5], d, n, m
	elseif number == 1
		n -= 1
		return [0.5, x1[2]], d, n, m
	elseif number == 3
		m -= 1
		return [x1[1], 0.5], d, n, m
	elseif error("Don't know which wall was hit")
	end
	
end


# Calculates the coordinates of places and of each hit obstacle of the collisions given the time
function collisions_classical(x, v, r, tmax, prec::Integer=64)

	set_bigfloat_precision(prec)
	#x = big(x); v = big(v); r = big(r)
	x = [BigFloat("$(x[1])"), BigFloat("$(x[2])")]; v = [BigFloat("$(v[1])"), BigFloat("$(v[2])")]; r = BigFloat("$r")
	places = Array{BigFloat, 1}[]
	circles = Vector{BigInt}[]
	speeds = Array{BigFloat, 1}[] # Modification to collect speeds too; may be turned off
	times = BigFloat[]
	# Put the starting point into the places array
	push!(places, x)	
	push!(speeds, v)
	push!(times, 0)

	# Initial square (n, m)
	n = ifloor(x[1] + 0.5)
	m = ifloor(x[2] + 0.5)
	
	# Place the first initial position into square [-0.5, 0.5)^2
	x -= [n, m]
	
	if norm(x) < r
		error("The initial position cannot be inside an obstacle")
	end

	t = 0
	while t <= tmax
		# Will hit or miss? Check the condition for hitting
		
		vcrossx = abs(crossz(v, x))
		vr = norm(v)*r
		discr = vr^2 - vcrossx^2
		
		if vcrossx < vr && (-dot(v, x) - sqrt(discr)) >= 0 # To make sure it does not go backwards
			#println("hit")
			# Then reflect

			t1 = (-dot(v, x) - sqrt(discr))/norm(v)^2
			
			N0 = x + v*t1
			N = N0/norm(N0)
			
			# Velocity, place and time immediately after the collision
			v1 = v - 2dot(v, N)*N
			x1 = x + v*t1
			t += t1
			
			push!(places, x1 + [n, m])
			push!(circles, [n, m])
			push!(speeds, v1)
			push!(times, t)
			#print("{$(r1[1] + n), $(r1[2] + m)}, ")
			
			x, d, n, m = crossing(x1, v1, n, m)
						
			# The speed direction will stay the same
			v = v1
			
			# The time will increment once again from collision point to the wall
			t += d/norm(v)
			
		else # If it misses the ball
			#println("miss")
			v1 = v
			x1 = x
			
			# Now hit the wall
			x, d, n, m = crossing(x1, v1, n, m)
						
			# The speed direction will stay the same
			v = v1
			
			# The time will increment
			t += d/norm(v)
		end
	end
	return places, circles, speeds, times
end


# Calculates the coordinates of obstacles hit at the first collision
function first_collision_classical(x::Vector, v::Vector, r::Real, precision::Integer=64)

	set_bigfloat_precision(precision)
	circles = Vector[]

	# Initial square (n, m)
	n = ifloor(x[1] + 0.5)
	m = ifloor(x[2] + 0.5)
	
	# Place the first initial position into square [-0.5, 0.5)^2
	x -= [n, m]
	
	if norm(x) < r
		error("The initial position cannot be inside an obstacle")
	end


	while abs(crossz(v, x)) > norm(v)*r # misses the ball

		v1 = v
		x1 = x
		
		# Now hit the wall
		x, d, n, m = crossing(x1, v1, n, m)
					
		# The speed direction will stay the same
		v = v1
	
	end
	
	(n, m)
end


### 3D collisions

function crossing3d(x, v, n, m, l)
	# Outputs new position in cube [-0.5, 0.5)^3 and updates the numbers of new square
	
	# Minimum positive times to nearest crossing - concise formulas p. 48 - 49
	array_times = BigFloat[]
	#k = v[2]/v[1]	
	for i = 1:3
		for j = 1:2
			push!(array_times, ((-1)^j*0.5 - x[i])/v[i])
		end	
	end
	
	# Extract the minimum positive time value
	minpos = Inf
	number = 0
	for i = 1:length(array_times)
		if array_times[i] < minpos && array_times[i] > 0
			minpos = array_times[i]
			number = i
		end				
	end
	
	t = array_times[number]
	x1 = x + v*t
	d = norm(v*t)
	
	# Possible outcomes depending on which wall will be hit
	if number == 2
		n += 1
		return [-0.5, x1[2], x1[3]], d, n, m, l
	elseif number == 4
		m += 1
		return [x1[1], -0.5, x1[3]], d, n, m, l
	elseif number == 1
		n -= 1
		return [0.5, x1[2], x1[3]], d, n, m, l
	elseif number == 3
		m -= 1
		return [x1[1], 0.5, x1[3]], d, n, m, l
	elseif number == 6
		l += 1
		return [x1[1], x1[2], -0.5], d, n, m, l
	elseif number == 5
		l -= 1
		return [x1[1], x1[2], 0.5], d, n, m, l
	elseif error("Don't know which wall was hit")
	end
	
end


function collisions3d_classical(x0, v0, r, tmax, precision::Integer=64)

	set_bigfloat_precision(precision)
	x0 = big(x0); v0 = big(v0); r = big(r)
	#x0 = [BigFloat("$(x0[1])"), BigFloat("$(x0[2])"), BigFloat("$(x0[3])")]; v0 = [BigFloat("$(v0[1])"), BigFloat("$(v0[2])"), BigFloat("$(v0[3])")]; r = BigFloat("$r")
	# Normalize speed
	v0 /= norm(v0)
	places = Vector{BigFloat}[]
	circles = Vector{BigInt}[]
	speeds = Vector{BigFloat}[]
	times = BigFloat[]
	# Put the starting point into the places array
	push!(places, x0)	
	push!(speeds, v0)
	push!(times, 0)

	# Initial square (n, m)
	n = ifloor(x0[1] + 0.5)
	m = ifloor(x0[2] + 0.5)
	l = ifloor(x0[3] + 0.5)
	
	# Place the first initial position into cube [-0.5, 0.5)^3
	x0 -= [n, m, l]
	
	if norm(x0) < r
		error("The initial position cannot be inside an obstacle")
	end

	t = 0
	while t <= tmax
		# Will hit or miss? Check the condition for hitting - derivation pp. 83-84
		xcrossv = norm(cross(v0, x0))
		vr = norm(v0)*r
		discr = vr^2 - xcrossv^2
		if xcrossv < vr && (-dot(v0, x0) - sqrt(discr)) > 0
			#println("hit")
			# Then reflect
			
			
			# Throw away complex time values if any, but there shouldn't be
			#if discr >= 0
			t1 = (-dot(v0, x0) - sqrt(discr))/norm(v0)^2
			#end
			
			N0 = x0 + v0*t1
			N = N0/norm(N0)
			
			# Velocity, place and time immediately after the collision
			v1 = v0 - 2*dot(v0, N)*N
			x1 = x0 + v0*t1
			t += t1
			
			push!(places, x1 + [n, m, l])
			push!(circles, [n, m, l])
			push!(speeds, v1)
			push!(times, t)
			norm(x1), norm(x1 + [n,m,l] - [n,m,l])
			#print("{$(r1[1] + n), $(r1[2] + m)}, ")
			
			#x0 += v0*0.001 # Don't know why it was here. Commented out.
			
			x0, d, n, m, l = crossing3d(x1, v1, n, m, l)
			
						
			# The speed direction will stay the same
			v0 = v1
			
			# The time will increment once again from collision point to the wall
			t += d/norm(v0)

			
		else # If it misses the ball
			#println("miss")
			v1 = v0
			x1 = x0
			
			# Now hit the wall
			x0, d, n, m, l = crossing3d(x1, v1, n, m, l)
						
			# The speed direction will stay the same
			v0 = v1
			
			# The time will increment
			t += d/norm(v0)
		end
	end
	return places, circles, speeds, times
end




function first_collision3d_classical(x0, v0, r, precision::Integer=64)

	set_bigfloat_precision(precision)
	x0 = big(x0); v0 = big(v0); r = big(r)
	#x0 = [BigFloat("$(x0[1])"), BigFloat("$(x0[2])"), BigFloat("$(x0[3])")]; v0 = [BigFloat("$(v0[1])"), BigFloat("$(v0[2])"), BigFloat("$(v0[3])")]; r = BigFloat("$r")
	# Normalize speed
	v0 /= norm(v0)

	# Initial square (n, m)
	n = ifloor(x0[1] + 0.5)
	m = ifloor(x0[2] + 0.5)
	l = ifloor(x0[3] + 0.5)
	
	# Place the first initial position into cube [-0.5, 0.5)^3
	x0 -= [n, m, l]
	
	if norm(x0) < r
		error("The initial position cannot be inside an obstacle")
	end

	t = 0

	# Will hit or miss? Check the condition for hitting - derivation pp. 83-84
	#xcrossv = norm(cross(v0, x0))
	#vr = norm(v0)*r
	#discr = vr^2 - xcrossv^2
	#if xcrossv < vr && (-dot(v0, x0) - sqrt(discr)) > 0
	while norm(cross(v0, x0)) > norm(v0)*r

		v1 = v0
		x1 = x0
		
		# Now hit the wall
		x0, d, n, m, l = crossing3d(x1, v1, n, m, l)
					
		# The speed direction will stay the same
		v0 = v1
		
		# The time will increment
		t += d/norm(v0)
	end
	
	return (n, m, l)
end


# End of module
end
