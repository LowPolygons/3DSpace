module coordsGen
	implicit none

contains
	subroutine generateFiles(success)
		implicit none

		integer :: numparticles, c1

		integer, dimension(8) :: seed
		integer :: nParticles
		real, dimension(3) :: axisDimension
		real, dimension(3) :: minBoundary
		real :: cutoff
		logical :: exists
		real, dimension(:,:), allocatable :: Aparticles
		real, dimension(:,:), allocatable :: particles

		logical :: success

		print *, "Number of particles: "
		read(*,*) nParticles

		print *, "Seed: "
		read(*,*) c1

		print *, "Cutoff: "
		read(*,*) cutoff

		print *, "Input minimum coordinate (enter one coordinate at a time and press enter): "
		read(*,*) minBoundary

		print *, "Input maximum coordinate (enter one coordinate at a time and press enter): "
		read(*,*) axisDimension

		axisDimension = axisDimension - minBoundary
		seed = c1

		allocate(Aparticles(nParticles,3))
		allocate(particles(nParticles,3))

		call random_seed(put=seed)
		call random_number(Aparticles)
		
		do c1 = 1, nParticles
			particles(c1, 1:3) = [&
				& minBoundary(1) + Aparticles(c1, 1)*axisDimension(1), &
				& minBoundary(2) + Aparticles(c1, 2)*axisDimension(2), &
				& minBoundary(3) + Aparticles(c1, 3)*axisDimension(3) &
			&]
		end do

		inquire(file="coordinates.txt", exist=exists)

		if (exists) then
			open(10, file="coordinates.txt", status="old")
		else
			open(10, file="coordinates.txt", status="new")
		end if

		9 format(f20.17, 4x, f20.17, 4x, f20.17)

		do c1 = 1, nParticles
			write(10, 9) particles(c1,1), particles(c1,2), particles(c1,3)
		end do

		close(10)

		inquire(file="config.txt", exist=exists)

		if (exists) then
			open(11, file="config.txt", status="old")
		else
			open(11, file="config.txt", status="new")
		end if 

		!NPARTICLES
		write(11, *) nParticles
		!CUTOFF
		write(11, *) cutoff
		!LOWERBOUND
		write(11, 9) minBoundary
		!UPPERBOUND
		write(11, 9) (minBoundary+axisDimension)

		close(11)

		print *, "Generated"
		print *, "........................."
		print *, "........................."
		deallocate(Aparticles, particles)

		inquire(file="coordinates.txt", exist=exists)
		success = .false.

		if (exists) then
			inquire(file="config.txt", exist=exists)
			if (exists) then
				success = .true.
			end if 
		end if 
	end subroutine generateFiles
	
end module coordsGen