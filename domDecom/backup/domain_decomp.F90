!gfortran -O2 -march=native vectors.F90 -o o2test
!test ./o2test

! A more intense implementation of MPI where the main space is divided among processors
! When running, use a power of 2 processors

!any non obvious vector functions go here

module coordsGen
	implicit none

contains
	function generateFiles() result(success)
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
		character :: logging
		logical :: success

		print *, "-> Number of particles: "
		read(*,*) nParticles

		print *, "-> Seed: "
		read(*,*) c1

		print *, "-> Cutoff: "
		read(*,*) cutoff

		print *, "-> Input minimum coordinate (enter one coordinate at a time and press enter): "
		read(*,*) minBoundary

		print *, "-> Input maximum coordinate (enter one coordinate at a time and press enter): "
		read(*,*) axisDimension

		print *, "-> Logging (Y/N, Case Sensitive): "
		read(*,*) logging

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
				& minBoundary(3) + Aparticles(c1, 3)*axisDimension(3)  &
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
		!LOGGING
		write(11, *) logging

		close(11)

		print "(a18//)", "-> Generated Files"
		deallocate(Aparticles, particles)

		inquire(file="coordinates.txt", exist=exists)
		success = .false.

		if (exists) then
			inquire(file="config.txt", exist=exists)
			if (exists) then
				success = .true.
			end if 
		end if 
	end function generateFiles
	
end module coordsGen

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
	use iso_fortran_env
	use coordsGen

	implicit none

	real, dimension(1:3) :: lowerBound, upperBound, currParticle
	real, dimension(6) :: curr2
	integer :: counter1, counter2, counter3
	integer :: i, j
	real :: cutoff, tempReal

	real, dimension(:,:), allocatable :: particlePositions

	integer(kind=int64) :: particlePairCount
	integer :: totalSum

	!mpi stuff
	integer :: ierr, rank, nprocs, numParticles, numRecieved, numRecievedAfter
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio, arraySize, oppositeCount, leftoverTracker
	integer, dimension(1:3) :: numDomains
	real, dimension(1:6) :: lowerAndUpper !boundary
	integer, dimension(1:6) :: adjacents
	real, dimension(:,:), allocatable :: processorRange, sentProcessData, processDataHolder !the 2nd : is going to be 6, lowerboundxyz, upperboundxyz
	real, dimension(:,:), allocatable :: transferThese, receivedThese, transferTheseAfter, receivedTheseAfter ,allDataReceived  !for use in each process to sort the data into the left/right, down/up, back/forward regions
	real, dimension(:,:), allocatable :: addAtEnd !for use with adjcanet process particle things
	!integer, dimension(:) , allocatable :: processTracker
	integer :: processTracker
	logical :: success, dologging
	character :: generateNewFiles, logging


	! lowerBound = [-0.10000000149011612, -18.150000001490117 , -14.540000001490116]![0.0, 0.0, 0.0] !
	! upperBound = [18.150000001490117, 0.10000000149011612, 0.10000000149011612]   ![1.0, 1.0, 1.0] !
	! cutoff = 4.0 !0.25 !4.0 
	
	totalSum = 0


	call MPI_INIT(ierr)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	


	if (rank == 0) then
		print *, "Generate new particles and config file (Y/N, Case Sensitive):"
		read (*,*) generateNewFiles

		success = .false.

		if (generateNewFiles == "Y") then
			success = generateFiles()
		else
			inquire(file="coordinates.txt", exist=success)

			if (success) then
				inquire(file="config.txt", exist=success)
			end if 
		end if

		if (.not. success) stop "-> Files not successfully created. Program ending"
	end if

	call MPI_BARRIER(MPI_COMM_WORLD, ierr)

	open(10, file="config.txt", status="old")
	
	read(10, *) numParticles
	read(10, *) cutoff
	read(10, *) lowerBound
	read(10, *) upperBound
	read(10, *) logging

	close(10)


	if (rank == 0) then
		dologging = .false.

		if (logging == "Y") dologging = .true.

		print *, ""
		print *, "<==> Program Config Parameters <==>"
		print *, "Cutoff:", cutoff
		print *, "Number of Particles:", numParticles
		print *, "Lower Bound:", lowerBound
		print *, "Upper Bound:", upperBound
		print *, "Logging Enabled:", dologging
		print *, "<======>"
		print "(//a24/)", "<==> Program Start <==>"

		if (dologging .and. rank == 0) print *, "--> Reading file data"

		!open(11, file="pData.dat", status="old")
		!open(11, file="particles_data.dat", status="old")
		open(11, file="coordinates.txt", status="old")
		!open(11, file="copperblock1.dat", status="old")

		allocate(particlePositions(numParticles, 3))

		do counter1 = 1, size(particlepositions)/3
			read(11, *) particlepositions(counter1,1), particlepositions(counter1,2), particlepositions(counter1,3)
		end do
		close(11)

		! Code that is good for vectorised manipulation of the particles
		! do counter1 = 1, size(particlepositions)/3
		! 	particlepositions(counter1,1) = particlepositions(counter1,1) - 1
		! 	particlepositions(counter1,2) = particlepositions(counter1,2) - 1
		! 	particlepositions(counter1,3) = particlepositions(counter1,3) - 1
		! end do

		allocate(processorRange(nprocs,6))

		if (dologging .and. rank == 0) print *, "--> File data read"

		numDomains = determineNoSplits(nprocs)

		cutoff = min( min(cutoff, (upperBound(1)-lowerBound(1))/ (numDomains(1)*2)), &
			min( (upperBound(2)-lowerBound(2))/(numDomains(2)*2), (upperBound(3)-lowerBound(3))/(numDomains(3)*2)) )

		12 format(" --> Space divided across:", i4, " processors, Number of domains per axis:", i4, ",", i4, "," i4)
		if (dologging .and. rank == 0) print 12, nprocs, numDomains(1), numDomains(2), numDomains(3)

		processorRange = getBounds(nprocs, lowerBound, upperBound, numDomains) !determineBounds(nprocs, lowerBound, upperBound) !

		allocate(processDataHolder(size(particlePositions)/3, 3))
		processDataHolder = upperBound(1)*100 !arbitrary number so we know when its wrong
		
		if (dologging .and. rank == 0) print *, "--> Distributing particles among processors"
		
		do counter2 = 1, nprocs !currentProcessr
			processTracker = 0
			processDataHolder = upperBound(1)*10
			curr2(1:6) = processorRange(counter2, 1:6)

			do counter1 = 1, size(particlePositions)/3 !currentParticle

				if (particlePositions(counter1, 1) <= upperBound(1)) then
					currParticle = particlePositions(counter1, 1:3)
					if (        currParticle(1) >= curr2(1) .and. currParticle(1) < curr2(4)&
						& .and. currParticle(2) >= curr2(2) .and. currParticle(2) < curr2(5)&
						& .and. currParticle(3) >= curr2(3) .and. currParticle(3) < curr2(6)) then
						processTracker = processTracker + 1
						processDataHolder(processTracker, 1:3) = currParticle(1:3)
						!tag the current particle so it can be glossed over
						particlePositions(counter1, 1:3) = upperBound(1)*1000
					end if
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

				!Particles
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


	!send to all other processors the num of domains
	call MPI_BCAST(numDomains, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(cutoff, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
	adjacents = getAdjacents(rank, numDomains)
	!adjacents = getAdjacantRegions(nprocs, rank)
	!print *, rank, adjacents

	if (dologging .and. rank == 0) print *, "--> Distribution Complete"
	

	particlePairCount = 0

	if (dologging .and. rank == 0) print *, "--> Comparing particles in own process domain."	
	do counter2 = 1, arraySize
		do counter3 = counter2+1, arraySize
			!  particlePairCount = particlePairCount + STSinRange(sentProcessData(counter2, 1:3),&
			!  &sentProcessData(counter3, 1:3), cutoff, rank)
			particlePairCount = particlePairCount + inRange(sentProcessData(counter2, 1:3),&
			& sentProcessData(counter3, 1:3), upperBound, lowerBound, cutoff)
		end do
	end do

	call MPI_BARRIER(MPI_COMM_WORLD, ierr)

	if (dologging .and. rank == 0) print *, "--> Own particle comparisons finished."	

	allocate(allDataReceived(arraySize, 3))
	allDataReceived = sentProcessData

	counter2 = 0

	! print *, rank, ":", lowerAndUpper
	if (dologging .and. rank == 0) print *, "--> Starting iterations."	
	do counter1 = 1, 5, 2

		counter2 = counter2 + 1 

		allocate(transferThese(size(allDataReceived)/3, 3))
		allocate(transferTheseAfter(size(allDataReceived)/3, 3)) 
		allocate(addAtEnd(size(allDataReceived)/3, 3))

		processTracker = 0 
		oppositeCount = 0 
		leftoverTracker = 0 


		do counter3 = 1, size(allDataReceived)/3
			currParticle = allDataReceived(counter3, 1:3)
			if ( currParticle(counter2) < lowerAndUpper(counter2)+cutoff ) then
				processTracker = processTracker + 1
				transferThese(processTracker, 1:3) = currParticle(1:3)
			else if ( currParticle(counter2) >= lowerAndUpper(counter2+3)-cutoff) then
				oppositeCount = oppositeCount + 1
				transferTheseAfter(oppositeCount, 1:3) = currParticle(1:3)
			else
				leftoverTracker = leftoverTracker + 1
				addAtEnd(leftoverTracker, 1:3) = currParticle(1:3)
			end if
		end do

		numRecieved = 0
		CALL MPI_Sendrecv(processTracker, 1, MPI_INTEGER, adjacents(counter1), 2, &
           				     numRecieved, 1, MPI_INTEGER, adjacents(counter1+1), 2, &
           		 MPI_COMM_WORLD, status1, ierr)

		numRecievedAfter = 0

		CALL MPI_Sendrecv(oppositeCount, 1, MPI_INTEGER, adjacents(counter1+1), 3, &
           			   numRecievedAfter, 1, MPI_INTEGER, adjacents(counter1), 3, &
            		   MPI_COMM_WORLD, status1, ierr)

		allocate(receivedThese(numRecieved, 3))
		allocate(receivedTheseAfter(numRecievedAfter, 3))

		! do counter3 = 1, processTracker
		! 	transferThese(counter3, 1:3) = transferThese(counter3, 1:3) - lowerAndUpper(1:3)
		! end do 

		! if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		! if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		! if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		! do counter3 = 1, oppositeCount
		! 	transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) - currParticle
		! end do 

		!sendreceive the particles to the lower and from the higher
		CALL MPI_Sendrecv(transferThese(1:processTracker, 1:3), processTracker*3, MPI_REAL, adjacents(counter1), 10, &
          					receivedThese(1:numRecieved, 1:3), numRecieved*3, MPI_REAL, adjacents(counter1+1), 10, &
          			MPI_COMM_WORLD, status1, ierr)

		!For the after particles
		call MPI_Sendrecv(transferTheseAfter(1:oppositeCount, 1:3), oppositeCount*3, MPI_REAL,  adjacents(counter1+1), 10, &
				    receivedTheseAfter(1:numRecievedAfter, 1:3), numRecievedAfter*3, MPI_REAL,  adjacents(counter1), 10, &
				& MPI_COMM_WORLD, status1, ierr)

		! if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		! if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		! if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		! do counter3 = 1, numRecieved
		! 	receivedThese(counter3, 1:3) = receivedThese(counter3, 1:3) + currParticle
		! end do

		! do counter3 = 1, numRecievedAfter
		! 	receivedTheseAfter(counter3, 1:3) = receivedTheseAfter(counter3, 1:3) + lowerAndUpper(1:3)
		! end do

		! do counter3 = 1, processTracker
		! 	transferThese(counter3, 1:3) = transferThese(counter3, 1:3) + lowerAndUpper(1:3)
		! end do 

		! if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		! if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		! if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		! do counter3 = 1, oppositeCount
		! 	transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) + currParticle
		! end do 

		do i = 1, size(sentProcessData)/3 
			do j = 1, size(receivedThese)/3
				!  particlePairCount = particlePairCount + STSinRange(sentProcessData(i, 1:3),&
				!  & receivedThese(j, 1:3), cutoff, rank)				
				particlePairCount = particlePairCount + inRange(sentProcessData(i, 1:3),&
				& receivedThese(j, 1:3), upperBound, lowerBound, cutoff)
			end do
		end do

		deallocate(allDataReceived)

		allocate(allDataReceived(processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter, 3))

		if (processTracker /= 0) allDataReceived(1:processTracker, 1:3) = transferThese(1:processTracker, 1:3)

		if (oppositeCount /= 0) allDataReceived(processTracker+1 : processTracker+oppositeCount, 1:3) &
		& = transferTheseAfter(1:oppositeCount, 1:3)

		if (leftoverTracker /= 0) allDataReceived(processTracker+oppositeCount+1 : processTracker+oppositeCount+leftoverTracker, 1:3) &
		& = addAtEnd(1:leftoverTracker, 1:3)

		if (numRecieved /= 0) allDataReceived(&
			&processTracker+oppositeCount+leftoverTracker+1: &
			& (processTracker+oppositeCount+leftoverTracker+numRecieved), &
			1:3) = receivedThese(1:numRecieved, 1:3)

		if (numRecievedAfter /= 0) allDataReceived(&
			&processTracker+oppositeCount+leftoverTracker+numRecieved+1: &
			& (processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter), &
			1:3) = receivedTheseAfter(1:numRecievedAfter, 1:3)


		deallocate(transferThese) 
		deallocate(transferTheseAfter) 
		deallocate(receivedThese) 
		deallocate(receivedTheseAfter) 
		deallocate(addAtEnd)

		43 format(" --> Iteration ", i1, " complete.")
		if (dologging .and. rank == 0) print 43, counter2
	end do 

	deallocate(allDataReceived)

	if (rank /= 0) then
		call MPI_SEND(particlePairCount, 1, MPI_INTEGER8, 0, 3, MPI_COMM_WORLD, ierr)
	else
		do i = 1, nprocs-1
			counter1 = 0
			call MPI_RECV(counter1, 1, MPI_INTEGER8, i, 3, MPI_COMM_WORLD, status1, ierr)
			particlePairCount = particlePairCount + counter1 
		end do
		print *, "<==> Total number of unique pairs with cutoff ", cutoff, " : ", particlePairCount, "<==>"
	end if


	call MPI_FINALIZE(ierr)

contains
	integer function STSinRange(p1, p2, cutoff, rank) result(count)
		!no periodic boundary
		real, dimension(1:3), intent(in) :: p1, p2
		real, dimension(1:3) :: temp
		real, intent(in):: cutoff
		integer, intent(in) :: rank

		temp = p2 - p1

		count = 0

		!if (vectorMod(temp) /= 0) print *, p1, p2, vectorMod(temp)
		if (vectorMod(temp) < cutoff) count = 1
	end function STSinRange

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

	function getAdjacents(currProcess, axisSplits) result(adjacents)
		integer, intent(in) :: currProcess
		integer, dimension(1:3) :: axisSplits, currPos, currPosDup, remainders, adders, periodicAdditions
		integer :: processDup, axis
		integer, dimension(6) :: adjacents
	
		currPos = 0
		processDup = currProcess + 1 !makes the beneath code work

		adders(1) = 1
		adders(2) = axisSplits(1)
		adders(3) = axisSplits(1)*axisSplits(2)
		!This next block of code determines which layer per axis the current processor is in

		!Start with Z - for clarification this is which Z layer it lies in, starting @ 1 not 0
		remainders(3) = processDup !for later
		do while(processDup > 0)
			currPos(3) = currPos(3) + 1
			processDup = processDup - adders(3)
		end do 
		!undoing the last subtrat which took it out of the loop
		processDup = processDup + adders(3)

		remainders(2) = processDup ! for later
		!DO the same for Y
		do while(processDup > 0)
			currPos(2) = currPos(2) + 1
			processDup = processDup - adders(2)
		end do 
		processDup = processDup + adders(2)

		remainders(1) = processDup
		!remainer is the x
		currPos(1) = processDup

		periodicAdditions(1) = remainders(2) + remainders(3)
		periodicAdditions(2) = remainders(1) + remainders(3)
		periodicAdditions(3) = remainders(1) + remainders(2)

		currPosDup = currPos !want to know the originals but also modify them 

		do axis = 1, 3
			if (axisSplits(axis) == 1) then
				adjacents((axis-1)*2 + 1: (axis-1)*2 + 2) = ((currPosDup(1)-1)*adders(1)) + &
				& ((currPosDup(2)-1)*adders(2)) + ((currPosDup(3)-1)*adders(3))
			else
				!lower
				if (currPos(axis)-2 < 0) then
					currPosDup(axis) = axisSplits(axis)
				else
					currPosDup(axis) = currPos(axis)-1 !it does another -1 later
				end if

				adjacents( (axis-1)*2 + 1) = ((currPosDup(1)-1)*adders(1)) + & 
				& ((currPosDup(2)-1)*adders(2)) + ((currPosDup(3)-1)*adders(3))

				!upper
				if (currPos(axis)+1 > axisSplits(axis)) then
					currPosDup(axis) = 1
				else
					currPosDup(axis) = currPos(axis)+1
				end if

				adjacents( (axis-1)*2 + 2) = ((currPosDup(1)-1)*adders(1)) + & 
				& ((currPosDup(2)-1)*adders(2)) + ((currPosDup(3)-1)*adders(3))
				
				currPosDup = currPos
			end if
		end do 
		!print *, "-->", currProcess, " has its adjacent processors."
	end function getAdjacents

	function getBounds(nprocs, lowerBound, upperBound, axisSplits) result(boundaryRanges)
		integer, intent(in) :: nprocs
		real, dimension(1:3), intent(in) :: lowerBound, upperBound
		real, dimension(nprocs,6) :: boundaryRanges
		real, dimension(1:3) :: difference, rangePerAxis
		integer, dimension(1:3), intent(in) :: axisSplits
		integer :: x, y, z, numIterations

		difference = upperBound - lowerBound
		!axisSplits = determineNoSplits(nprocs)
		
		rangePerAxis(1) = difference(1) / axisSplits(1)
		rangePerAxis(2) = difference(2) / axisSplits(2)
		rangePerAxis(3) = difference(3) / axisSplits(3)

		!print *, "range per axis: ", rangePerAxis

		x = 0
		y = 0
		z = 0
		!Now loop through all possible regions in binary order and assign to appropriate boundaryRanges place
		do numIterations = 1, nprocs
			if (x >= axisSplits(1)) then
				y = y + 1
				x = 0
			end if
			if (y >= axisSplits(2)) then
				z = z + 1
				x = 0
				y = 0
			end if 
			!this should count up in binary
			boundaryRanges(numIterations, 1:3) = lowerBound + [x*rangePerAxis(1), y*rangePerAxis(2), z*rangePerAxis(3)]
			boundaryRanges(numIterations, 4:6) = lowerBound + [(x+1)*rangePerAxis(1), (y+1)*rangePerAxis(2), (z+1)*rangePerAxis(3)]

			x = x + 1
		end do
		print *, "--> Get Bounds Complete"
	end function getBounds 

	function determineNoSplits(nprocs) result(numDomains)
		integer, intent(in) :: nprocs
		integer :: i, currFactor, numCopy
		integer, dimension(3) :: numDomains
		integer :: numFactors
		integer, dimension(nprocs) :: factors

		factors = 0
		numFactors = 0
		numCopy = nprocs
		currFactor = 1

		!breaks a given number into a product of factors (not necessarily prime, often is though)
		do while (numCopy /= 1)
			currFactor = currFactor + 1
			!print *, currFactor, numCopy
			if ( mod(numCopy, currFactor) == 0) then
				numFactors = numFactors + 1
				factors(numFactors) = currFactor

				numCopy = numCopy / currFactor
				currFactor = 1
				!print *, currFactor, numCopy
			end if
		end do

		numFactors = numFactors + 1
		factors(numFactors) = 1

		numDomains = 1

		!split the factors among the domains
		do i = 0, numFactors-1
			numDomains(1+mod(i, 3)) = numDomains(1+mod(i, 3)) * factors(i+1)
		end do 

		!print *, nprocs, ":", numDomains

		print *, "--> Determine Number Of Splits done"	
	end function determineNoSplits

	real (kind=8) function timestamp() result(elapsed) 
		implicit none
		real(kind=8) :: start_time=0.0, end_time = 0.0 
		logical :: initialised = .false.

		if (initialised) then
			start_time = end_time 
			end_time = MPI_WTIME()
			elapsed = end_time - start_time
		else
			start_time = MPI_WTIME()
			initialised = .true.
		elapsed = 0.0
		end if
	end function timestamp

end program vectors