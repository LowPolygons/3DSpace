!gfortran -O2 -march=native vectors.F90 -o o2test
!test ./o2test

!any non obvious vector functions go here
module vectorFunctions
	implicit none
contains
	function vectorMod(v1) result(ret_mod)
		implicit none
		real, dimension(1:3), intent(in) :: v1
		real :: ret_mod
		
		ret_mod = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
	end function vectorMod
	
end module vectorFunctions

program vectors
	use vectorFunctions
	
	implicit none

	real, dimension(1:3) :: lowerBound, upperBound
	integer :: counter1, counter2, counter3, nParts
	real :: cutoff

	real, dimension(:,:), allocatable :: particlePositions

	integer, dimension(:), allocatable :: particlePairCount
	integer :: totalSum

	lowerBound = [0.0, 0.0, 0.0] ![-0.10000000149011612, -18.150000001490117 , -14.540000001490116]!
	upperBound = [1.0, 1.0, 1.0] ![18.150000001490117, 0.10000000149011612, 0.10000000149011612]   !
	cutoff = 0.25 !4.0 
	totalSum = 0

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Receiving position data from dat file
	!!! And populating position/velocity arrays 

	!Positions
	!Not strictly necessary
	9 format(f20.17, 4x, f20.17, 4x, f20.17)

	open(11, file="pData.dat", status="old")
	!open(11, file="copperblock1.dat", status="old")

	read(11, *) nParts! particle size

	allocate(particlePositions(nParts,3))
	allocate(particlePairCount(nParts))


	do counter1 = 1, size(particlePositions)/3
		particlePairCount(counter1) = 0
	end do

	do counter1 = 1, size(particlepositions)/3
		read(11, 9) particlepositions(counter1,1), particlepositions(counter1,2), particlepositions(counter1,3)
	end do

	do counter1 = 1, size(particlePositions)/3
		particlepositions(counter1,1) = particlepositions(counter1,1)
		particlepositions(counter1,2) = particlepositions(counter1,2)
		particlepositions(counter1,3) = particlepositions(counter1,3)
	end do

	close(11)

	! Loops through time steps in the case of velocity involvement
	! Calculates all the pairs of particles that need checking and runs the inRange function
	do counter2 = 1, size(particlePositions)/3
		do counter3 = counter2, size(particlePositions)/3
			!shouldn't check itself
			if (counter3 /= counter2) then
				!print *, counter2, counter3
				particlePairCount(counter2) = particlePairCount(counter2) + inRange(particlePositions(counter2, 1:3),&
				&particlePositions(counter3, 1:3), upperBound, lowerBound, cutoff)
				!print *, "||||||||||||||||||||||||||||||||||"
			end if
		end do
	end do

	do counter1 = 1, size(particlePositions)/3
		!print *, counter1, particlePairCount(counter1)
		totalSum = totalSum + particlePairCount(counter1)
	end do

	print *, "Total number of unique interactions: ", totalSum

	deallocate(particlePositions)
	deallocate(particlePairCount)

contains
	!from the current and target particle, it checks the quickest 3 direct routes in each axis and checks the mod of that vector
	function inRange(pos1, pos2, upperBound, lowerBound, cutoff) result(success)
		implicit none
		real, dimension(1:3), intent(in) :: pos1, pos2
		real, dimension(1:3) :: upperBound, lowerBound
		real, intent(in) :: cutoff
		integer :: success
		real, dimension(1:3) :: temp

		success = 0
		!compare the distance between the two, the distance from p1Left and p2Right, and p1Right and p2Left
		!these match up the only 3 direct routes from each particle to the other and so the shortest of these and then Modulused
		!will tell you what is in cut off
		temp(1) = min( abs(pos2(1)-pos1(1)),&
		 &min( (	(pos1(1)-lowerBound(1))	+ (upperBound(1)-pos2(1))), &
		 &(	(pos2(1)-lowerBound(1))	+ (upperBound(1)-pos1(1))) ) )

		temp(2) = min( abs(pos2(2)-pos1(2)),&
		 & min( (	(pos1(2)-lowerBound(2))	+ (upperBound(2)-pos2(2))), &
		 &(	(pos2(2)-lowerBound(2))	+ (upperBound(2)-pos1(2))) ) )
		
		temp(3) = min( abs(pos2(3)-pos1(3)),&
		 &min( (	(pos1(3)-lowerBound(3))	+ (upperBound(3)-pos2(3))), &
		 &(	(pos2(3)-lowerBound(3))	+ (upperBound(3)-pos1(3))) ) )

		if (vectorMod(temp) < cutoff) success = 1
	end function inRange


end program vectors