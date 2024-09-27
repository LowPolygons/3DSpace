!gfortran -O2 -march=native vectors.F90 -o o2test
!test ./o2test

! A more intense implementation of MPI where the main space is divided among processors
! When running, use a power of 2 processors

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

	real, dimension(1:3) :: lowerBound, upperBound, currParticle
	real, dimension(6) :: curr2
	integer :: counter1, counter2, counter3
	integer :: i, j
	real :: cutoff, tempReal
	integer :: request

	real, dimension(:,:), allocatable :: particlePositions

	integer :: particlePairCount
	integer :: totalSum, tempVar

	!mpi stuff
	integer :: ierr, rank, nprocs, numParticles, numRecieved, numRecievedAfter
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio, arraySize, oppositeCount, leftoverTracker

	real, dimension(1:6) :: lowerAndUpper !boundary
	integer, dimension(1:6) :: adjacents
	real, dimension(:,:), allocatable :: processorRange, sentProcessData, processDataHolder !the 2nd : is going to be 6, lowerboundxyz, upperboundxyz
	real, dimension(:,:), allocatable :: transferThese, receivedThese, transferTheseAfter, receivedTheseAfter ,allDataReceived  !for use in each process to sort the data into the left/right, down/up, back/forward regions
	real, dimension(:,:), allocatable :: adParticles, addAtEnd !for use with adjcanet process particle things
	!integer, dimension(:) , allocatable :: processTracker
	integer :: processTracker
	lowerBound = [0.0, 0.0, 0.0]  !
	upperBound = [1.0, 1.0, 1.0] !
	cutoff = 0.25
	totalSum = 0
	tempVar = 0

	! do counter1 = 1, size(particlePositions)/3
	! 	particlePairCount(counter1) = 0
	! end do
	call MPI_INIT(ierr)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	


	if (rank == 0) then
		4 format ("[", f20.17, ",",f20.17, ",",f20.17, " ]")
		!test to see data is loaded in. it is
		!do counter1 = 1, size(particlePositions)/3
		!	print 4, particlePositions(counter1,1),particlePositions(counter1,2), particlePositions(counter1,3)
		!end do

		!send only the desired data to each process. Only open in msater process
		open(11, file="pData.dat", status="old")

		read(11, *) numParticles

		allocate(particlePositions(numParticles, 3))

		do counter1 = 1, size(particlepositions)/3
			read(11, *) particlepositions(counter1,1), particlepositions(counter1,2), particlepositions(counter1,3)
		end do

		do counter1 = 1, size(particlePositions)/3
			particlepositions(counter1,1) = particlepositions(counter1,1)
			particlepositions(counter1,2) = particlepositions(counter1,2)
			particlepositions(counter1,3) = particlepositions(counter1,3)
		end do

		close(11)
		

		allocate(processorRange(nprocs,6))

		processorRange = determineBounds(nprocs, lowerBound, upperBound)
		!the range of each process
		allocate(processDataHolder(size(particlePositions)/3, 3))
		!allocate(processTracker(nprocs))
		processDataHolder = upperBound(1)*100 !arbitrary number so we know when its wrong
	
		!finding which particles go where and sending this to the various process'
		do counter2 = 1, nprocs !currentProcessr
			processTracker = 0
			processDataHolder = upperBound(1)*10
			curr2(1:6) = processorRange(counter2, 1:6)

			do counter1 = 1, size(particlePositions)/3 !currentParticle
				currParticle = particlePositions(counter1, 1:3)
				
				if (       currParticle(1) >= curr2(1) .and. currParticle(1) < curr2(4)&
					&.and. currParticle(2) >= curr2(2) .and. currParticle(2) < curr2(5)&
					&.and. currParticle(3) >= curr2(3) .and. currParticle(3) < curr2(6)) then
					processTracker = processTracker + 1
					processDataHolder(processTracker, 1:3) = currParticle(1:3)
				end if
			end do 

			!now you know the number of particles, its boundary, and its particle data, send them
			if (counter2 == 1) then
				arraySize = processTracker
				lowerAndUpper(1:6) = curr2(1:6)
				allocate(sentProcessData(arraySize, 3))

				sentProcessData(1:arraySize, 1:3) = processDataHolder(1:arraySize, 1:3)
			else
				!numparticles
				call MPI_SEND(processTracker, 1, MPI_INTEGER, counter2-1, counter2-1, MPI_COMM_WORLD, ierr)
				!boundary
				call MPI_SEND(curr2(1:6), 6, MPI_REAL, counter2-1, (counter2-1)*2, MPI_COMM_WORLD, ierr)
				!numParticles
				call MPI_SEND(processDataHolder(1:processTracker, 1:3), 3*processTracker, MPI_REAL, counter2-1, &
				&(counter2-1)*3, MPI_COMM_WORLD, ierr)
			end if
		end do
	end if 

	!send the array size and bounds to each process
	if (rank /= 0) then
		call MPI_RECV(arraySize, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, status1, ierr)
		call MPI_RECV(lowerAndUpper, 6, MPI_REAL, 0, rank*2, MPI_COMM_WORLD, status1, ierr)

		allocate(sentProcessData(arraySize, 3))

		call MPI_RECV(sentProcessData, 3*arraySize, MPI_REAL, 0, rank*3, MPI_COMM_WORLD, status1, ierr)
	end if
	
	!determiens what processors are left and right from the current
	adjacents = getAdjacantRegions(nprocs, rank)

	print *, rank, arraySize

	!local pair counting
	particlePairCount = 0
	do counter2 = 1, arraySize
		do counter3 = counter2+1, arraySize
			particlePairCount = particlePairCount + inRange(sentProcessData(counter2, 1:3),&
			&sentProcessData(counter3, 1:3), cutoff, rank)
		end do
	end do

	call MPI_BARRIER(MPI_COMM_WORLD, ierr)

	!this array is going to be used to store all relevant particles so that previosuly received ones trasnfer across loops,
	allocate(allDataReceived(arraySize, 3))
	allDataReceived = sentProcessData

	counter2 = 0
	!if (rank == 2) print *, lowerAndUpper
	!this will represent how many particles you are sending to the adjacent process
	do counter1 = 1, 5, 2!2 !1, 3, 5 ; these are the adjacents which you send data to, add one for the one you receive from

		counter2 = counter2 + 1 !represents which axis you are on
		!first, establish the particles you are sending

		allocate(transferThese(size(allDataReceived)/3, 3)) !temporary storage
		allocate(transferTheseAfter(size(allDataReceived)/3, 3)) !temporary storage
		allocate(addAtEnd(size(allDataReceived)/3, 3))

		processTracker = 0 !for transferThese
		oppositeCount = 0 !for transferTheseAfter
		leftoverTracker = 0 !for leftovers

		do counter3 = 1, size(allDataReceived)/3
			currParticle = allDataReceived(counter3, 1:3)

			!if it lies in the lower half of the region accessible by the lower bound + cutoff, it might be reachable in adjacent
			if ( currParticle(counter2) < lowerAndUpper(counter2)+cutoff ) then
				processTracker = processTracker + 1
				transferThese(processTracker, 1:3) = currParticle(1:3)
			else if ( currParticle(counter2) > lowerAndUpper(counter2+3)-cutoff) then
				oppositeCount = oppositeCount + 1
				transferTheseAfter(oppositeCount, 1:3) = currParticle(1:3)
			else
				leftoverTracker = leftoverTracker + 1
				addAtEnd(leftoverTracker, 1:3) = currParticle(1:3)
			end if
		end do


		call MPI_ISEND(processTracker, 1, MPI_INTEGER, adjacents(counter1), rank*1*counter2, MPI_COMM_WORLD, request, ierr)
		numRecieved = 0
		!determine how many particles you're about to receive
		call MPI_RECV(numRecieved, 1, MPI_INTEGER, adjacents(counter1+1), adjacents(counter1+1)*1*counter2, &
		& MPI_COMM_WORLD, status1, ierr)

		!now wait to ensure it arrived
		call MPI_WAIT(request, status1, ierr)

		!now the second
		call MPI_ISEND(oppositeCount, 1, MPI_INTEGER, adjacents(counter1), rank*5*counter2, MPI_COMM_WORLD, request, ierr)
		numRecievedAfter = 0

		call MPI_RECV(numRecievedAfter, 1, MPI_INTEGER, adjacents(counter1+1), adjacents(counter1+1)*5*counter2, &
		&MPI_COMM_WORLD, status1, ierr)

		!now wait to ensure it arrived
		call MPI_WAIT(request, status1, ierr)

		!now ready to receive the particles
		allocate(receivedThese(numRecieved, 3))
		allocate(receivedTheseAfter(numRecievedAfter, 3))

		!now, ensure the particles are the offsets rather than their world space positions to allow for periodic boundaries to function
		do counter3 = 1, processTracker
			!Sending left/down/backward therefore find the distance from the lowest coordinate
			transferThese(counter3, 1:3) = transferThese(counter3, 1:3) - lowerAndUpper(1:3)
			!if (rank == 0) print*, counter2, transferThese(counter3, 1:3)
		end do 

		!for after particles
		!Reuse currentParticle as the offset you are subtracting from each particle
		if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		do counter3 = 1, oppositeCount
			!Sending left/down/backward therefore find the distance from the lowest coordinate
			transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) - currParticle
		end do 

		!send the actual particles
		call MPI_ISEND(transferThese(1:processTracker, 1:3), processTracker*3, MPI_REAL, &
		& adjacents(counter1), rank*2*counter2, MPI_COMM_WORLD, request, ierr)

		!Receive the particles
		call MPI_RECV(receivedThese(1:numRecieved, 1:3), numRecieved*3, MPI_REAL, &
		& adjacents(counter1+1), adjacents(counter1+1)*2*counter2, MPI_COMM_WORLD, status1, ierr)

		!now wait to ensure it arrived
		call MPI_WAIT(request, status1, ierr)

		!For the after particles
		call MPI_ISEND(transferTheseAfter(1:oppositeCount, 1:3), oppositeCount*3, MPI_REAL, &
		& adjacents(counter1), rank*7*counter2, MPI_COMM_WORLD, request, ierr)

		call MPI_RECV(receivedTheseAfter(1:numRecievedAfter, 1:3), numRecievedAfter*3, MPI_REAL, &
		& adjacents(counter1+1), adjacents(counter1+1)*7*counter2, MPI_COMM_WORLD, status1, ierr)


		!now wait to ensure it arrived
		call MPI_WAIT(request, status1, ierr)

		!Particles are now in the correct coordinate space. From here, add this to the vertex of the bound which has the upper value in the current axis

		!Reuse currentParticle as the offset you are adding to each particle
		if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		do counter3 = 1, numRecieved
			receivedThese(counter3, 1:3) = receivedThese(counter3, 1:3) + currParticle
		end do

		!Particles after need to just be added to the bottom left corner always
		do counter3 = 1, numRecievedAfter
			receivedTheseAfter(counter3, 1:3) = receivedTheseAfter(counter3, 1:3) + lowerAndUpper(1:3)
		end do

		!undoing the calculation before sending it so that the original coordinates get saved
		do counter3 = 1, processTracker
			!Sending left/down/backward therefore find the distance from the lowest coordinate
			transferThese(counter3, 1:3) = transferThese(counter3, 1:3) + lowerAndUpper(1:3)
			!if (rank == 0) print*, counter2, transferThese(counter3, 1:3)
		end do 

		!undoing transferTheseAfter
		if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		do counter3 = 1, oppositeCount
			!Sending left/down/backward therefore find the distance from the lowest coordinate
			transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) + currParticle
		end do 

		!Now the particles are in the correct space, you do the pair searching

		do i = 1, arraySize !only want to loop through seld to self particles
			do j = 1, numRecieved
				particlePairCount = particlePairCount + inRange(sentProcessData(i, 1:3),&
				& receivedThese(j, 1:3), cutoff, rank)
			end do
		end do

		!Lastly, receive the particle data from the opposite side as this ensures there are no missed pairs

		!transfering data into allDataReceived for furtuer loops
		deallocate(allDataReceived)
		allocate(allDataReceived(processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter, 3))
		!allocate(allDataReceived(arraySize+numRecieved, 3))

		!allDataReceived(1:arraySize, 1:3) = sentProcessData(1:arraySize, 1:3)
		!add transferThese, transferAfter, addAtEnd here
		if (processTracker > 0) allDataReceived(1:processTracker, 1:3) = transferThese(1:processTracker, 1:3)

		if (oppositeCount > 0) allDataReceived(processTracker+1 : processTracker+oppositeCount, 1:3) &
		& = transferTheseAfter(1:oppositeCount, 1:3)

		if (leftoverTracker > 0) allDataReceived(processTracker+oppositeCount+1 : processTracker+oppositeCount+leftoverTracker, 1:3) &
		& = addAtEnd(1:leftoverTracker, 1:3)

		if (numRecieved > 0) allDataReceived(&
			&processTracker+oppositeCount+leftoverTracker+1: &
			& (processTracker+oppositeCount+leftoverTracker+numRecieved), &
			1:3) = receivedThese(1:numRecieved, 1:3)

		if (numRecievedAfter > 0) allDataReceived(&
			&processTracker+oppositeCount+leftoverTracker+numRecieved+1: &
			& (processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter), &
			1:3) = receivedTheseAfter(1:numRecievedAfter, 1:3)

		! allDataReceived((arraySize+numRecieved)+1:(arraySize+numRecieved+numRecievedAfter), 1:3)&
		! & = receivedTheseAfter(1:numRecievedAfter, 1:3)

		deallocate(transferThese) 
		deallocate(transferTheseAfter) 
		deallocate(receivedThese) 
		deallocate(receivedTheseAfter) 
		deallocate(addAtEnd)
	end do 

	!print *, rank, particlePairCount
	deallocate(allDataReceived)

	if (rank /= 0) then
		call MPI_SEND(particlePairCount, 1, MPI_INT, 0, rank*20, MPI_COMM_WORLD, ierr)
	else
		do i = 1, nprocs-1
			counter1 = 0
			call MPI_RECV(counter1, 1, MPI_INT, i, i*20, MPI_COMM_WORLD, status1, ierr)
			particlePairCount = particlePairCount + counter1 
		end do
		print *, "Total number of Particles: ", particlePairCount
	end if


	call MPI_FINALIZE(ierr)

contains
	integer function inRange(p1, p2, cutoff, rank) result(count)
		!no periodic boundary
		real, dimension(1:3), intent(in) :: p1, p2
		real, dimension(1:3) :: temp
		real, intent(in):: cutoff
		integer, intent(in) :: rank

		temp = p2 - p1

		count = 0
		if (vectorMod(temp) < cutoff) count = 1
	end function inRange

	!its important to note that each process
	function getAdjacantRegions(nprocs, currRank) result(adjacents)
		integer, dimension(6) :: adjacents !xL xR, yU yD, zF zB
		integer, intent(in) :: nprocs, currRank
		integer :: totalPowerOfTwo, curr, i
		integer :: numXSlices, numYSlices, numZSlices, usedProcs
		integer, dimension(3) :: slicesToCurr, duplicateSlices

		totalPowerOfTwo = 0 !2^0 is 1 so if only 1 processor is put in, it does it all

		do i = 1, nprocs !nprocs is fine because one scales linearly and one scales exponentially
			if (2**i > nprocs) then
				exit
			else
				totalPowerOfTwo = i
			end if
		end do 

		numXSlices = 0
		numYSlices = 0
		numZSlices = 0

		!determines how many times it has been cut in each direction
		usedProcs = 2**totalPowerOfTwo - 1 !as its 0 to 7 for eg not 1 to 8
		do i = 0, usedProcs, 4
			if (i+4 <= usedProcs) numZSlices = numZSlices + 1
		end do 
		usedProcs = usedProcs - numZSlices*4

		do i = 0, usedProcs , 2
			if (i+2 <= usedProcs) numYSlices = numYSlices + 1
		end do 
		usedProcs = usedProcs - numYSlices*2
		
		numXSlices = usedProcs
		!as above

		!now do the same for the curr Process
		slicesToCurr = 0
		curr = currRank
		do i = 0, curr, 4
			if (i+4 <= curr) slicesToCurr(3) = slicesToCurr(3) + 1
		end do 
		curr = curr - slicesToCurr(3)*4

		do i = 0, curr, 2
			if (i+2 <= curr)slicesToCurr(2) = slicesToCurr(2) + 1
		end do 
		curr = curr - slicesToCurr(2)*2


		slicesToCurr(1) = curr

		!if (currRank == 6) print *, slicesToCurr
		!now have the number of slices in in each direction to get to your process
		!just going to do it line by line initially cus thats the easiest way
		duplicateSlices = slicesToCurr

		!left
		duplicateSlices(1) = duplicateSlices(1) - 1
		if (duplicateSlices(1) < 0) duplicateSlices(1) = numXSlices

		adjacents(1) = (duplicateSlices(1) + duplicateSlices(2)*2 + duplicateSlices(3)*4)
		duplicateSlices(1) = slicesToCurr(1)

		!right
		duplicateSlices(1) = duplicateSlices(1) + 1
		if (duplicateSlices(1) > numXSlices) duplicateSlices(1) = 0

		adjacents(2) = (duplicateSlices(1) + duplicateSlices(2)*2 + duplicateSlices(3)*4)
		duplicateSlices(1) = slicesToCurr(1)

		!down
		duplicateSlices(2) = duplicateSlices(2) - 1
		if (duplicateSlices(2) < 0) duplicateSlices(2) = numYSlices

		adjacents(3) = (duplicateSlices(1) + duplicateSlices(2)*2 + duplicateSlices(3)*4)
		duplicateSlices(2) = slicesToCurr(2)

		!up
		duplicateSlices(2) = duplicateSlices(2) + 1
		if (duplicateSlices(2) > numYSlices) duplicateSlices(2) = 0

		adjacents(4) = (duplicateSlices(1) + duplicateSlices(2)*2 + duplicateSlices(3)*4)
		duplicateSlices(2) = slicesToCurr(2)

		!backward
		duplicateSlices(3) = duplicateSlices(3) - 1
		if (duplicateSlices(3) < 0) duplicateSlices(3) = numZSlices

		adjacents(5) = (duplicateSlices(1) + duplicateSlices(2)*2 + duplicateSlices(3)*4)
		duplicateSlices(3) = slicesToCurr(3)

		!forward
		duplicateSlices(3) = duplicateSlices(3) + 1
		if (duplicateSlices(3) > numZSlices) duplicateSlices(3) = 0

		adjacents(6) = (duplicateSlices(1) + duplicateSlices(2)*2 + duplicateSlices(3)*4)
		duplicateSlices(3) = slicesToCurr(3)

	end function getAdjacantRegions

	function determineBounds(nprocs, lowerBound, upperBound) result(boundaryRanges)
		!this will return an array of 6 wide vectors for each processor, outlining what their boundaries shall be
		!this program will be designed to only use a power of 2 number of processors and it will disregard any others
		!it will do logical divisions, halving the total X, then Y then Z in a cycle for each power of two the number is
		integer, intent(in) :: nprocs
		real, dimension(3), intent(in) :: lowerBound, upperBound
		real, dimension(3) :: maxBoundRange
		integer, dimension(3) :: divsPerAxis !by default this will be the difference in lower vs upper bound,
		!divsPerAxis will say how many times each axis have been cut, for easier looping later
		integer :: totalPowerOfTwo, i, j, k, numProcessorsUsed, numAssignments
		real, dimension(nprocs,6) :: boundaryRanges !this will store the ranges of processors, in order of Rank 0 -> max

		totalPowerOfTwo = 0 !2^0 is 1 so if only 1 processor is put in, it does it all

		do i = 1, nprocs !nprocs is fine because one scales linearly and one scales exponentially
			if (2**i > nprocs) then
				exit
			else
				totalPowerOfTwo = i
			end if
		end do 

		do i = 1, nprocs
			boundaryRanges(i, 1:6) = [0.0,0.0,0.0,0.0,0.0,0.0]
		end do

		numProcessorsUsed = 2**totalPowerOfTwo !this is how many of the process' will be used. Not a superb long term 
		!solution but it makes the divisions easy

		!now, it will begin the process of determine each boxes' max bounds range
		maxBoundRange = upperBound - lowerBound
		!the number of slices needed is the total powers of two

		divsPerAxis = 1 !ignore if equal to zero !not zero because we want the loop to run at least once and if it doesnt get sliced the loop wont run

		do i = 1, totalPowerOfTwo
			j = max(1, mod(i, totalPowerOfTwo+1)) !so it repeats and is 1 and totalPowerOfTwo inclusive
			! tested in lua in lua demo website curr = math.max(1, (curr+1) % (maxCount + 1))
			!print *, i, j
			maxBoundRange(j) = maxBoundRange(j) / 2
			divsPerAxis(j) = divsPerAxis(j) + 1
		end do 


		numAssignments = 0
		i = 0
		j = 0
		k = 0

		!based on the previous calcualted values, it determines what the range of values each process will make use of
		do while(numAssignments < numProcessorsUsed)
			i = max(1, mod(i+1, divsPerAxis(1)+1))
			if (i == 1) j = max(1, mod(j+1, divsPerAxis(2)+1))
			if (i == 1 .and. j == 1) k = max(1, mod(k+1, divsPerAxis(3)+1))
			!assigning values
			numAssignments = numAssignments + 1


			boundaryRanges(numAssignments, 1) = maxBoundRange(1)*(i-1)
			boundaryRanges(numAssignments, 2) = maxBoundRange(2)*(j-1)
			boundaryRanges(numAssignments, 3) = maxBoundRange(3)*(k-1)
			boundaryRanges(numAssignments, 4) = maxBoundRange(1)*(i)
			boundaryRanges(numAssignments, 5) = maxBoundRange(2)*(j)
			boundaryRanges(numAssignments, 6) = maxBoundRange(3)*(k)
		end do

		! 100 format("Lower : [", 3f10.7, "], Upper : [" 3f10.7, "]")
		! do i = 1, nprocs
		! 	print 100, boundaryRanges(i,1:3),boundaryRanges(i,4:6)
		! end do
		!!!!!!!!!!!!!!!!!!!!result
		! mpirun -n 8 ./test
		! Lower : [ 0.0000000 0.0000000 0.0000000], Upper : [ 0.5000000 0.5000000 0.5000000]
		! Lower : [ 0.5000000 0.0000000 0.0000000], Upper : [ 1.0000000 0.5000000 0.5000000]
		! Lower : [ 0.0000000 0.5000000 0.0000000], Upper : [ 0.5000000 1.0000000 0.5000000]
		! Lower : [ 0.5000000 0.5000000 0.0000000], Upper : [ 1.0000000 1.0000000 0.5000000]
		! Lower : [ 0.0000000 0.0000000 0.5000000], Upper : [ 0.5000000 0.5000000 1.0000000]
		! Lower : [ 0.5000000 0.0000000 0.5000000], Upper : [ 1.0000000 0.5000000 1.0000000]
		! Lower : [ 0.0000000 0.5000000 0.5000000], Upper : [ 0.5000000 1.0000000 1.0000000]
		! Lower : [ 0.5000000 0.5000000 0.5000000], Upper : [ 1.0000000 1.0000000 1.0000000]
	end function determineBounds

end program vectors