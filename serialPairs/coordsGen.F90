program main
	implicit none

	integer :: numparticles, c1
	real, dimension(3) :: range
	real, dimension(10000,3) :: particles, velocities

	call random_number(particles)

	open(10, file="pData.dat", status="new")

	1 format(f20.17, 4x, f20.17, 4x, f20.17)

	do c1 = 1, 10000
		write(10, 1) particles(c1,1), particles(c1,2), particles(c1,3)
	end do

	close(10)

	call random_number(velocities)

	open(11, file="vData.dat", status="new")

	do c1 = 1, 10000
		write(11, 1) velocities(c1,1), velocities(c1,2), velocities(c1,3)
	end do

	close(11)

	! open(11, file="pData.dat", status="old")

	! do c1 = 1, 100
	! 	read(11, 1) newParticles(c1,1), newParticles(c1,2), newParticles(c1,3)
	! end do

	! close(11)

	! print *, "Old: ", particles
	! print *, "New: ", newParticles

end program main