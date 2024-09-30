program main
	implicit none

	integer :: numparticles, c1
	real, dimension(3) :: range
	real, dimension(10000,3) :: Aparticles
	integer, dimension(10000,3) :: particles

	call random_number(Aparticles)
	
	particles(1:10000, 1:3) = Aparticles(1:10000, 1:3)*2

	open(10, file="IData.dat", status="new")


	do c1 = 1, 10000
		write(10, *) particles(c1,1), particles(c1,2), particles(c1,3)
	end do

	close(10)


	! open(11, file="pData.dat", status="old")

	! do c1 = 1, 100
	! 	read(11, 1) newParticles(c1,1), newParticles(c1,2), newParticles(c1,3)
	! end do

	! close(11)

	! print *, "Old: ", particles
	! print *, "New: ", newParticles

end program main