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

!in the event of particles moving around, the boundary function can ensure it enters periodic space properl
module boundaries
	use vectorFunctions
	implicit none
contains
	function checkBoundary(currPosAfterSim, lowerBound, upperBound) result(newPos)
		implicit none
		real, dimension(1:3), intent(in) :: currPosAfterSim, lowerBound, upperBound
		real, dimension(1:3) :: newPos
		real :: tempVariable
		!xyz handled separately as its smarter to do so	
		newPos = currPosAfterSim
		!x	
		if (currPosAfterSim(1) < lowerBound(1)) then
			tempVariable = currPosAfterSim(1) - lowerBound(1)
			newPos(1) = upperBound(1) + tempVariable
		else if (currPosAfterSim(1) > upperBound(1)) then
			tempVariable = currPosAfterSim(1) - upperBound(1)
			newPos(1) = lowerBound(1) + tempVariable		
		end if	
		!y
		if (currPosAfterSim(2) < lowerBound(2)) then
			tempVariable = currPosAfterSim(2) - lowerBound(2)
			newPos(2) = upperBound(2) + tempVariable
		else if (currPosAfterSim(2) > upperBound(2)) then
			tempVariable = currPosAfterSim(2) - upperBound(2)
			newPos(2) = lowerBound(2) + tempVariable		
		end if
		!z
		if (currPosAfterSim(3) < lowerBound(3)) then
			tempVariable = currPosAfterSim(2) - lowerBound(3)
			newPos(3) = upperBound(3) + tempVariable
		else if (currPosAfterSim(3) > upperBound(3)) then
			tempVariable = currPosAfterSim(3) - upperBound(3)
			newPos(3) = lowerBound(3) + tempVariable		
		end if		
	end function checkBoundary
end module boundaries


program vectors
	use vectorFunctions
	use boundaries
	use mpi
	
	implicit none

	real, dimension(1:3) :: lowerBound, upperBound
	integer :: numOfParticles
	integer :: counter1, counter2, counter3
	real :: cutoff, randoVar
	integer :: pairInteractionsThisTick

	real, dimension(:), allocatable :: particlePositionsTotalX, particlePosPermX, particlePosDynX
	real, dimension(:), allocatable :: particlePositionsTotalY, particlePosPermY, particlePosDynY
	real, dimension(:), allocatable :: particlePositionsTotalZ, particlePosPermZ, particlePosDynZ

	integer, dimension(:), allocatable :: particlePairCount
	integer :: totalSum, tempVar

	!mpi stuff
	integer :: ierr, rank, nprocs
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio, arraySize
	
	lowerBound = [0.0, 0.0, 0.0]  
	upperBound =  [1.0, 1.0, 1.0] 
	cutoff = 0.5
	totalSum = 0
	tempVar = 0

	!initlaise mpi 
	call MPI_INIT(ierr)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	

	!initialise the numberOfparticles
	if (rank == 0) then
		open(11, file="pData.dat", status="old")

		read(11, *) numOfParticles

		allocate(particlePositionsTotalX(numOfParticles))
		allocate(particlePositionsTotalY(numOfParticles))
		allocate(particlePositionsTotalZ(numOfParticles))
		! allocate(particlePosPerm(numOfParticles,3))
		! allocate( particlePosDyn(numOfParticles,3))

		do counter1 = 1, numOfParticles
			read(11, *) particlePositionsTotalX(counter1), particlePositionsTotalY(counter1), particlePositionsTotalZ(counter1)
		end do 

		arraySize = (numOfParticles / nprocs) + mod(numOfParticles, nprocs)

	end if

	!each processor now knows how large to make their array
	call MPI_BCAST(arraySize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

	allocate(particlePosPermX(arraySize))
	allocate( particlePosDynX(arraySize))
	allocate(particlePosPermY(arraySize))
	allocate( particlePosDynY(arraySize))
	allocate(particlePosPermZ(arraySize))
	allocate( particlePosDynZ(arraySize))

	particlePosPermX = 5
	 particlePosDynX = 5
	particlePosPermY = 5
	 particlePosDynY = 5
	particlePosPermZ = 5
	 particlePosDynZ = 5 

		
	!sending particle data
	if (rank == 0) then
		!handles Rank0 separately
		counter2 = (numOfParticles / nprocs) !default integer division 
		counter3 = mod(numOfParticles, nprocs)
		
		particlePosPermX(1:counter2) = particlePositionsTotalX(1:counter2)
		particlePosPermY(1:counter2) = particlePositionsTotalY(1:counter2)
		particlePosPermZ(1:counter2) = particlePositionsTotalZ(1:counter2)

		do counter1 = 1, nprocs-1
			if (counter1 == nprocs-1) then
				!particlePositionsTotal(counter1*counter2 + 1:(counter1+1)*counter2 + counter3, 3)

				call MPI_SEND(particlePositionsTotalX(counter1*counter2:(counter1+1)*counter2 + counter3), &
				&3*(counter2+counter3), MPI_REAL, counter1, rank, MPI_COMM_WORLD, ierr) !*3 cus of the 2D elemnt of it
			else
				!particlePositionsTotal(counter1*counter2 + 1:(counter1+1)*counter2, 3)
				call MPI_SEND(particlePositionsTotalX(counter1*counter2 + 1:(counter1+1)*counter2), 3*counter2, &
				&MPI_REAL, counter1, rank, MPI_COMM_WORLD, ierr)
			end if
		end do
	else 
		if (rank == nprocs-1) then
			call MPI_RECV(particlePosPermX(1:arraySize), arraySize*3, MPI_REAL, 0, 0, MPI_COMM_WORLD, status1, ierr)
		else
			call MPI_RECV(particlePosPermX(1:arraySize), arraySize*3, MPI_REAL, 0, 0, MPI_COMM_WORLD, status1, ierr)
		end if 
	end if

	print *, rank, size(particlePosPermX)/3, particlePosPermX(1)

	call MPI_FINALIZE(ierr)
contains
	!self explanatory
	function updateParticlePos(pos, vel) result(ret_pos)
		implicit none
		real, dimension(1:3) ::  ret_pos
		real, dimension(1:3), intent(in) :: pos, vel
		
		ret_pos = pos + vel
	end function updateParticlePos

	!from the current and target particle, it checks all regions immediately adjacent to the main space and checks if any of the particles are in range
	function inRange(pos1, pos2, upperBound, lowerBound, cutoff) result(success)
		implicit none
		real, dimension(1:3), intent(in) :: pos1, pos2
		real, dimension(1:3) :: upperBound, lowerBound
		real, intent(in) :: cutoff
		integer :: success
		real, dimension(1:3) :: temp, distanceFromLower, disFromLow 
		logical, dimension(3,2) :: activeSpaces
		logical, dimension(-1:1,-1:1,-1:1) :: spaceMultipliers
		integer :: cX, cY, cZ

		disFromlow = pos1 - lowerBound
		distanceFromLower = pos2 - lowerBound
		spaceMultipliers = .true.
		!six operations to try optimise
		!checking if the cut off distance even reaches the neighbouring regions 

		if ( (lowerBound(1) + -1*(upperBound(1)-lowerBound(1)) + distanceFromLower(1) + cutOff) < pos1(1) ) then
			spaceMultipliers(-1 ,-1:1, -1:1) = .false.
		end if

		 if ( (lowerBound(1) + 1*(upperBound(1)-lowerBound(1)) + distanceFromLower(1) - cutOff) > pos1(1) ) then
			spaceMultipliers(1 ,-1:1, -1:1) = .false.
		 end if

		 if ( (lowerBound(2) + -1*(upperBound(2)-lowerBound(2)) + distanceFromLower(2) + cutOff) < pos1(2) ) then
		 	spaceMultipliers(-1:1 ,-1, -1:1) = .false.
		 end if

		 if ( (lowerBound(2) + 1*(upperBound(2)-lowerBound(2)) + distanceFromLower(2) - cutOff) > pos1(2) ) then
		 	spaceMultipliers(-1:1, 1, -1:1) = .false.
		 end if

		 if ( (lowerBound(3) + -1*(upperBound(3)-lowerBound(3)) + distanceFromLower(3) + cutOff) < pos1(3) ) then
		 	spaceMultipliers(-1:1 ,-1:1, -1) = .false.
		 end if

		 if ( (lowerBound(3) + 1*(upperBound(3)-lowerBound(3)) + distanceFromLower(3) - cutOff) > pos1(3) ) then
		 	spaceMultipliers(-1:1, -1:1, 1) = .false.
		 end if

		success = 0
		!if (CUTOFF GREATER THAN 0.5)
		do cX = -1, 1, 1
			do cY = -1, 1, 1
				do cZ = -1, 1, 1
					if (spaceMultipliers(cX,cY,cZ)) then
						temp(1) = lowerBound(1) + cX*(upperBound(1)-lowerBound(1))
						temp(2) = lowerBound(2) + cY*(upperBound(2)-lowerBound(2))
						temp(3) = lowerBound(3) + cZ*(upperBound(3)-lowerBound(3))
						if (vectorMod((temp+distanceFromLower)-pos1) < cutoff) then
							success = success + 1
						end if		
					end if		
				end do
			end do
		end do	
		!if (CUTOFF LESS THAN 0.5)
	end function inRange



end program vectors