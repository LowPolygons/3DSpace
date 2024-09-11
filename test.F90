

program vectors
	implicit none
	
	real, dimension(1:3) :: p1 = [1.0, 2.0, 3.0]
	real, dimension(1:3) :: p2 = [4.0, 5.0, 6.0]
	real, dimension(1:3) :: p3
	
	p3 = vectorAdd(p1, p2)
	
	print "(3f5.1)", p3

	
	p3 = vectorSub(p1, p2)
	
	print "(3f5.1)", p3	
	
	p3 = vectorMul(p1, p2)
	
	print "(3f5.1)", p3	
	
	p3 = vectorDiv(p1, p2)
	
	print "(3f5.1)", p3

contains
	function vectorAdd(v1, v2) result(ret_vector)
		implicit none
		real, dimension(1:3), intent(in) :: v1, v2
		real, dimension(1:3) :: ret_vector
		
		integer :: counter
		
		do counter = 1, 3
			ret_vector(counter) = v1(counter) + v2(counter)
		end do
		
	end function vectorAdd
	
	function vectorSub(v1, v2) result(ret_vector)
		implicit none
		real, dimension(1:3), intent(in) :: v1, v2
		real, dimension(1:3) :: ret_vector
		
		integer :: counter
		
		do counter = 1, 3
			ret_vector(counter) = v1(counter) - v2(counter)
		end do
		
	end function vectorSub

	function vectorMul(v1, v2) result(ret_vector)
		implicit none
		real, dimension(1:3), intent(in) :: v1, v2
		real, dimension(1:3) :: ret_vector
		
		integer :: counter
		
		do counter = 1, 3
			ret_vector(counter) = v1(counter) * v2(counter)
		end do
		
	end function vectorMul
	
	
	function vectorDiv(v1, v2) result(ret_vector)
		implicit none
		real, dimension(1:3), intent(in) :: v1, v2
		real, dimension(1:3) :: ret_vector
		
		integer :: counter
		
		do counter = 1, 3
			ret_vector(counter) = v1(counter) / v2(counter)
		end do
		
	end function vectorDiv
	
	
	
end program vectors


