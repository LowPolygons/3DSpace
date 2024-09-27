!This will not work if NumParticles % nprocs != 0

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
	integer :: numOfParticles
	integer :: counter1, counter2, counter3, currProcessorTarget, nextReceivingFrom
	real :: cutoff, randoVar
	integer :: pairInteractionsThisTick

	real, dimension(:,:), allocatable :: particlePositionsTotal, particlePosPerm, particlePosDyn

	integer :: particlePairCount, selfToSelfCounter
	integer :: totalSum, tempVar

	!mpi stuff
	integer :: ierr, rank, nprocs, allProcsDone
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio, arraySize
	integer :: STStotal, DBTTotal !self to self total, divide by two total
	
	lowerBound = [0.0, 0.0, 0.0]  
	upperBound =  [1.0, 1.0, 1.0] 
	cutoff = 0.5
	totalSum = 0
	tempVar = 0

	!initlaise mpi 
	call MPI_INIT(ierr)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	


! Initialises the particle data from the files and sends the data to each processor. Ensure nparticles%nprocs = 0
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
	!initialise the numberOfparticles
	if (rank == 0) then
		open(11, file="pData.dat", status="old")

		read(11, *) numOfParticles

		allocate(particlePositionsTotal(numOfParticles,3))
		! allocate(particlePosPerm(numOfParticles,3))
		! allocate( particlePosDyn(numOfParticles,3))

		do counter1 = 1, numOfParticles
			read(11, *) particlePositionsTotal(counter1, 1), particlePositionsTotal(counter1, 2), particlePositionsTotal(counter1, 3)
		end do 

		arraySize = (numOfParticles / nprocs) + mod(numOfParticles, nprocs)

	end if

	!each processor now knows how large to make their array
	call MPI_BCAST(arraySize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

	allocate(particlePosPerm(arraySize,3))
	allocate( particlePosDyn(arraySize,3))

	particlePosPerm = 5!Some arbirarily large number so we can isolate what isnt wanted later
	particlePosDyn = 5 !Some arbirarily large number so we can isolate what isnt wanted later

		
	!sending particle data
	if (rank == 0) then
		!handles Rank0 separately
		counter2 = (numOfParticles / nprocs) !default integer division 
		counter3 = mod(numOfParticles, nprocs)
		
		particlePosPerm(1:counter2, 1:3) = particlePositionsTotal(1:counter2, 1:3)
		do counter1 = 1, nprocs-1
			if (counter1 == nprocs-1) then
				!particlePositionsTotal(counter1*counter2 + 1:(counter1+1)*counter2 + counter3, 3)

				call MPI_SEND(particlePositionsTotal(counter1*counter2 + 1:(counter1+1)*counter2 + counter3, 1:3), &
				&3*(counter2+counter3), MPI_REAL, counter1, rank, MPI_COMM_WORLD, ierr) !*3 cus of the 2D elemnt of it
			else
				!particlePositionsTotal(counter1*counter2 + 1:(counter1+1)*counter2, 3)
				call MPI_SEND(particlePositionsTotal(counter1*counter2 + 1:(counter1+1)*counter2, 1:3), 3*counter2, &
				&MPI_REAL, counter1, rank, MPI_COMM_WORLD, ierr)
			end if
		end do
	else 
		if (rank == nprocs-1) then
			call MPI_RECV(particlePosPerm(1:arraySize, 1:3), arraySize*3, MPI_REAL, 0, 0, MPI_COMM_WORLD, status1, ierr)
		else
			call MPI_RECV(particlePosPerm(1:arraySize, 1:3), arraySize*3, MPI_REAL, 0, 0, MPI_COMM_WORLD, status1, ierr)
		end if 
	end if


	! print *, "Rank: ", rank, "Length: ", size(particlePosPerm)/3
	!Initialising the Dynamic/Changing Array for the first time
	particlePosDyn = particlePosPerm

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////!
	!Compare particles to itself
	currProcessorTarget = mod((rank + 1), (nprocs))
	nextReceivingFrom = rank - 1
	if (nextReceivingFrom < 0) nextReceivingFrom = nprocs-1
	!print *, nextReceivingFrom, rank, currProcessorTarget

	particlePairCount = 0
	selfToSelfCounter = 0
	do counter1 = 1, size(particlePosPerm)/3
		do counter2 = counter1+1, size(particlePosPerm)/3
			selfToSelfCounter = selfToSelfCounter + inRange(particlePosPerm(counter1, 1:3),&
			&particlePosDyn(counter2, 1:3), upperBound, lowerBound, cutoff)
		end do
	end do

	!after checking self to self, report back to master process, as these results only appear once but the rest appear twice so these shouldnt be divided.
	if (rank /= 0) then
		call MPI_SEND(selfToSelfCounter, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, ierr)
	else
		STStotal = selfToSelfCounter
		do counter1 = 1, nprocs-1
			call MPI_RECV(selfToSelfCounter, 1, MPI_INTEGER, counter1, counter1, MPI_COMM_WORLD, status1, ierr)
			STStotal = STStotal + selfToSelfCounter
		end do

		! print *, STStotal
	end if

	!TOTRY: Send a start signal to begin the pipeline process

	!Executing the Pipeline
	do counter1 = 1, nprocs-1 !the counter is not directly referencing the current processor
		!send to the next processor in the pipeline the current data in the Dyn
		call MPI_SEND(particlePosDyn, size(particlePosDyn), MPI_REAL, currProcessorTarget, rank, MPI_COMM_WORLD, ierr) !sending your current data to the next

		!then overwrite the Dyn with the data from the previous processor in the pipeline
		call MPI_RECV(particlePosDyn, size(particlePosDyn), MPI_REAL, &
		&nextReceivingFrom, nextReceivingFrom, MPI_COMM_WORLD, status1, ierr)

		!pair counting
		do counter2 = 1, size(particlePosPerm)/3
			do counter3 = 1, size(particlePosPerm)/3
				particlePairCount = particlePairCount + inRange(particlePosPerm(counter2, 1:3),&
				&particlePosDyn(counter3, 1:3), upperBound, lowerBound, cutoff)
			end do
		end do

		!TOTRY: Send a signal to master process to confirm that the current interation is done
		!Possibly the Send process' are getting blocked because one may send and then hang whilst another process is still doing some calculations
		! if (rank /= 0) then
		! 	call MPI_SEND(1, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, ierr)
		! else
		! 	allProcsDone = 0
		! 	counter3 = 0
		! 	do counter2 = 1, nprocs-1
		! 		call MPI_RECV(counter3, 1, MPI_INTEGER, counter2, counter2, MPI_COMM_WORLD, status1, ierr)
		! 		allProcsDone = allProcsDone + 1
		! 	end do
		! end if
	end do

	print *, Rank, " Done"


	!Each processor now sends it's total DBT particle count to the master processor 
	!This takes them all in, divides them by two, and adds it to the self-to-self to get the final result
	if (rank /= 0) then
		call MPI_SEND(particlePairCount, 1, MPI_INTEGER, 0, rank+nprocs+1, MPI_COMM_WORLD, ierr)
	else
		DBTTotal = particlePairCount
		do counter1 = 1, nprocs-1
			call MPI_RECV(particlePairCount, 1, MPI_INTEGER, counter1, counter1+nprocs+1, MPI_COMM_WORLD, status1, ierr)
			DBTTotal = DBTTotal + particlePairCount
		end do

		DBTTotal = DBTTotal / 2 !divide by two as name suggests


		print *, "Total Number of Pairs: ", (DBTTotal + STStotal)
	end if	

	call MPI_FINALIZE(ierr)
contains

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