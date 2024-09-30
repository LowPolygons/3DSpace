program messingWithVals
	implicit none
	
	integer, dimension(10,2) :: radius
	integer, dimension(10,2) :: newRadius
	
	integer :: i, j
	
	do i = 1, 10
		radius(i,1) = i
		radius(i,2) = i**2
	end do

	open(10, file="numbers.dat", status="old")

	do i = 1, 10
		write(10, "(2i6)") radius(i,1), radius(i,2)
	end do

	close(10)

	open(11, file="numbers.dat", status="old")
	
	do i = 1, 10
		read(11, "(2i6)") newRadius(i,1), newRadius(i,2)
	end do

	do i = 1, 10
		print "(2i6)", newRadius(i,1), newRadius(i,2)
	end do

end program messingWithVals


