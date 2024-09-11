program main
	implicit none

	integer :: numparticles, c1
	real, dimension(3) :: range
	real, dimension(100000,3) :: particles, newParticles

	call random_number(particles)

	open(10, file="pData.dat", status="old")

	1 format(f20.17, 4x, f20.17, 4x, f20.17)

	do c1 = 1, 100000
		write(10, 1) particles(c1,1), particles(c1,2), particles(c1,3)
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