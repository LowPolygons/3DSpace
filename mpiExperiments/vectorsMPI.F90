!gfortran -O2 -march=native vectors.F90 -o o2test
!test ./o2test

! A basic implementation of MPI in the particles system, where in a system of N particles and P processors, each processor does N/P of the particles

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
	use mpi
	
	implicit none

	real, dimension(1:3) :: lowerBound, upperBound
	integer :: timeSteps = 1
	integer :: counter1, counter2, counter3
	real :: cutoff, randoVar
	integer :: pairInteractionsThisTick

	real, dimension(1000,3) :: particlePositions

	integer, dimension(1000) :: particlePairCount
	integer :: totalSum, tempVar, sts

	!mpi stuff
	integer :: ierr, rank, nprocs
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio
	integer :: myLower, myUpper

	lowerBound = [0.0, 0.0, 0.0]  ![-0.10000000149011612, -18.150000001490117 , -14.540000001490116]!
	upperBound =  [1.0, 1.0, 1.0] ![18.150000001490117, 0.10000000149011612, 0.10000000149011612]   !
	cutoff = 0.05
	totalSum = 0
	tempVar = 0
	sts = 0

	do counter1 = 1, size(particlePositions)/3
		particlePairCount(counter1) = 0
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Receiving position data from dat file
	!!! And populating position/velocity arrays 

	!Positions
	!Not strictly necessary
	!9 format(f20.17, 4x, f20.17, 4x, f20.17)

	open(11, file="pData.dat", status="old")
	!open(11, file="copperblock1.dat", status="old")

	
	do counter1 = 1, size(particlepositions)/3
		read(11, *) particlepositions(counter1,1), particlepositions(counter1,2), particlepositions(counter1,3)
	end do

	! do counter1 = 1, size(particlePositions)/3
		! particlepositions(counter1,1) = particlepositions(counter1,1)
		! particlepositions(counter1,2) = particlepositions(counter1,2)
		! particlepositions(counter1,3) = particlepositions(counter1,3)
	! end do

	close(11)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!initlaise mpi 
	call MPI_INIT(ierr)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	


	!each processor does an equal amount.
	!integer divide the number of total particles by the number of processors
	!also mod divide the number of total particles
	!if the rank is zero, add to the number of operations the mod result so any processor count works
	particleProcRatio = (size(particlePositions) / 3 / nprocs) !already truncated
	
	!dont consider mods initally
	!if (rank == 0) particleProcRatio = particleProcRatio + mod(size(particlePositions) / 3, nprocs)

	!counter2 wants to loop through its section of the particle counts
	!if (rank == 0) print *, particleProcRatio

	myLower =  1 + (rank)*particleProcRatio
	myUpper =  (rank+1)*particleProcRatio

	if (rank == nprocs-1) myUpper = myUpper + mod(size(particlePositions) / 3, nprocs)
	print *, rank, myLower, myUpper

	!loop through STS first
	! do counter2 = myLower, myUpper
	! 	do counter3 = myLower, myUpper
	! 		!shouldn't check itself
	! 		if (counter3 /= counter2) then
	! 			sts = sts + inRange(particlePositions(counter1, 1:3),&
	! 			&particlePositions(counter2, 1:3), upperBound, lowerBound, cutoff)
	! 		end if
	! 	end do
	! end do

	!send sts back to master
	! if (rank /= 0) then
	! 	call MPI_SEND(sts, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, ierr)
	! else
	! 	do counter1 = 1, nprocs-1
	! 		counter2 = 0
	! 		call MPI_RECV(counter2, 1, MPI_INT, counter1, counter1, MPI_COMM_WORLD, status1, ierr)
	! 		sts = sts + counter2
	! 	end do
	! 	print *, "Self to self: ", sts
	! end if 

	!you divide these by two
	!checking the modulus remainder in case nprocs isnt a factor of th particle count
	do counter1 = myLower, myUpper
		do counter2 = 1, size(particlePositions)/3
			if (counter1 /= counter2) then
				totalSum = totalSum + inRange(particlePositions(counter1, 1:3),&
				&particlePositions(counter2, 1:3), upperBound, lowerBound, cutoff)
			end if
		end do 
	end do 
	
	! do counter1 = 1, size(particlePositions)/3
	! 	totalSum = totalSum + particlePairCount(counter1)
	! end do

	totalSum = totalSum/2

	if (rank /= 0) call MPI_SEND(totalSum, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, ierr)


	if (rank == 0) then 
		do counter1 = 1, nprocs-1
			tempVar = 0
			call MPI_RECV(tempVar, 1, MPI_INT, counter1, counter1, MPI_COMM_WORLD, status1, ierr)
			totalSum = totalSum + tempVar
		end do
		print *, "Total unique pairs: ", totalSum
	end if

	call MPI_FINALIZE(ierr)
contains
	!self explanatory
	function updateParticlePos(pos, vel) result(ret_pos)
		implicit none
		real, dimension(1:3) ::  ret_pos
		real, dimension(1:3), intent(in) :: pos, vel
		
		ret_pos = pos + vel
	end function updateParticlePos

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
	!from the current and target particle, it checks all regions immediately adjacent to the main space and checks if any of the particles are in range
	! function inRange(pos1, pos2, upperBound, lowerBound, cutoff) result(success)
		! implicit none
		! real, dimension(1:3), intent(in) :: pos1, pos2
		! real, dimension(1:3) :: upperBound, lowerBound
		! real, intent(in) :: cutoff
		! integer :: success
		! real, dimension(1:3) :: temp, distanceFromLower, disFromLow 
		! logical, dimension(3,2) :: activeSpaces
		! logical, dimension(-1:1,-1:1,-1:1) :: spaceMultipliers
		! integer :: cX, cY, cZ

		! disFromlow = pos1 - lowerBound
		! distanceFromLower = pos2 - lowerBound
		! spaceMultipliers = .true.
		!six operations to try optimise
		!checking if the cut off distance even reaches the neighbouring regions 

		! if ( (lowerBound(1) + -1*(upperBound(1)-lowerBound(1)) + distanceFromLower(1) + cutOff) < pos1(1) ) then
			! spaceMultipliers(-1 ,-1:1, -1:1) = .false.
		! end if

		 ! if ( (lowerBound(1) + 1*(upperBound(1)-lowerBound(1)) + distanceFromLower(1) - cutOff) > pos1(1) ) then
			! spaceMultipliers(1 ,-1:1, -1:1) = .false.
		 ! end if

		 ! if ( (lowerBound(2) + -1*(upperBound(2)-lowerBound(2)) + distanceFromLower(2) + cutOff) < pos1(2) ) then
		 	! spaceMultipliers(-1:1 ,-1, -1:1) = .false.
		 ! end if

		 ! if ( (lowerBound(2) + 1*(upperBound(2)-lowerBound(2)) + distanceFromLower(2) - cutOff) > pos1(2) ) then
		 	! spaceMultipliers(-1:1, 1, -1:1) = .false.
		 ! end if

		 ! if ( (lowerBound(3) + -1*(upperBound(3)-lowerBound(3)) + distanceFromLower(3) + cutOff) < pos1(3) ) then
		 	! spaceMultipliers(-1:1 ,-1:1, -1) = .false.
		 ! end if

		 ! if ( (lowerBound(3) + 1*(upperBound(3)-lowerBound(3)) + distanceFromLower(3) - cutOff) > pos1(3) ) then
		 	! spaceMultipliers(-1:1, -1:1, 1) = .false.
		 ! end if

		! success = 0
	!	if (CUTOFF GREATER THAN 0.5)
		! do cX = -1, 1, 1
			! do cY = -1, 1, 1
				! do cZ = -1, 1, 1
					! if (spaceMultipliers(cX,cY,cZ)) then
						! temp(1) = lowerBound(1) + cX*(upperBound(1)-lowerBound(1))
						! temp(2) = lowerBound(2) + cY*(upperBound(2)-lowerBound(2))
						! temp(3) = lowerBound(3) + cZ*(upperBound(3)-lowerBound(3))
						! if (vectorMod((temp+distanceFromLower)-pos1) < cutoff) then
							! success = success + 1
						! end if		
					! end if		
				! end do
			! end do
		! end do	
		!if (CUTOFF LESS THAN 0.5)
	! end function inRange



end program vectors