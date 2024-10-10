! Domain Decomposition
! Usage is very self explanatory

! Module for generating the files you want
module coordsGen
	implicit none

contains
	function generateFiles() result(success)
		implicit none

		integer :: numparticles, c1

		integer, dimension(8) :: seed
		integer :: nParticles
		double precision, dimension(3) :: axisDimension
		double precision, dimension(3) :: minBoundary
		double precision :: cutoff
		logical :: exists
		double precision, dimension(:,:), allocatable :: unit_particles
		double precision, dimension(:,:), allocatable :: particles
		character :: logging
		logical :: success

		! User inputs
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

		! The main program uses an upper and lower bound, however this wants a minimum and a dimension as it makes changing the regions better
		axisDimension = axisDimension - minBoundary
		seed = c1 			! Not directly stored into seed as it will need 8 inputs every time which is a headache and unecessary for this program

		allocate(unit_particles(nParticles,3))
		allocate(particles(nParticles,3))

		call random_seed(put=seed)
		call random_number(unit_particles)
		
		!Formatting the particles to be distributed evenly across the given upper and lower bound
		do c1 = 1, nParticles
			particles(c1, 1:3) = [&
				& minBoundary(1) + unit_particles(c1, 1)*axisDimension(1), &
				& minBoundary(2) + unit_particles(c1, 2)*axisDimension(2), &
				& minBoundary(3) + unit_particles(c1, 3)*axisDimension(3)  &
			&]
		end do

		!Writing the particles to the coordinates file
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

		! Writing the parameters to the config file: note, these are not labelled in any produced file
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
		deallocate(unit_particles, particles)


		! Confirm the success of writing to files
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


! Vector functions
module vectorFunctions
	implicit none
contains
	!Self explanatory
	function vectorMod(v1) result(ret_mod)
		implicit none
		double precision, dimension(1:3), intent(in) :: v1
		double precision :: ret_mod
		
		ret_mod = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
	end function vectorMod
	
end module vectorFunctions

program vectors
	use vectorFunctions
	use mpi
	use iso_fortran_env
	use coordsGen

	implicit none

	!!Output-Relevant, fixed sized arrays
	double precision, dimension(1:3) :: lowerBound, upperBound	
	double precision, dimension(6) :: curr2 		!RAM-functionality variables
	double precision, dimension(3) :: currParticle 	!RAM-functionality variables
	integer, dimension(1:3) :: numDomains
	double precision, dimension(1:6) :: lowerAndUpper
	integer, dimension(1:6) :: adjacents						!What processors are Left,Right,Down,Up,Forward,Backward

	!Counter variables
	integer :: counter1, counter2, counter3, i, j

	!Output-Relevant variables
	double precision :: cutoff
	integer(kind=int64) :: particlePairCount
	integer :: processTracker
	logical :: success, dologging
	character :: generateNewFiles, logging
	integer :: particleProcRatio, arraySize, oppositeCount, leftoverTracker

	!Output-Relevant, variable sized arrays
	double precision, dimension(:,:), allocatable :: particlePositions										!all particle positions
	double precision, dimension(:,:), allocatable :: processorRange, sentProcessData, processDataHolder
	double precision, dimension(:,:), allocatable :: transferThese, receivedThese
	double precision, dimension(:,:), allocatable :: transferTheseAfter, receivedTheseAfter, allDataReceived 
	double precision, dimension(:,:), allocatable :: addAtEnd 

	!mpi stuff
	integer :: ierr, rank, nprocs, numParticles, numRecieved, numRecievedAfter
    integer, dimension(MPI_STATUS_SIZE) :: status1

	!Init for MPI
	call MPI_INIT(ierr)

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	

	! Confirming if new data should be created or not. Checks if relevant files exist, if not program ends with message
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

	! Master-process only things
	if (rank == 0) then
		!Initialiser for logging
		dologging = .false.
		if (logging == "Y") dologging = .true.

		!Priting config parameters
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

		open(11, file="coordinates.txt", status="old")

		! All particle positions, only known for this process at this time
		allocate(particlePositions(numParticles, 3))

		do counter1 = 1, size(particlepositions)/3
			read(11, *) particlepositions(counter1,1), particlepositions(counter1,2), particlepositions(counter1,3)
		end do

		close(11)

		if (dologging .and. rank == 0) print *, "--> File data read"

		! Number of slices in each axis to distribute processors
		numDomains = determineNoSplits(nprocs)

		! In case an improper cutoff is put in, checks it
		cutoff = min( min(cutoff, (upperBound(1)-lowerBound(1))/ (numDomains(1)*2)), &
			min( (upperBound(2)-lowerBound(2))/(numDomains(2)*2), (upperBound(3)-lowerBound(3))/(numDomains(3)*2)) )

		12 format(" --> Space divided across:", i4, " processors, Number of domains per axis:", i4, ",", i4, "," i4)
		if (dologging .and. rank == 0) print 12, nprocs, numDomains(1), numDomains(2), numDomains(3)

		!The upper and lower bound of each process
		allocate(processorRange(nprocs,6))
		processorRange = getBounds(nprocs, lowerBound, upperBound, numDomains)

		!This will store the particles for the process that is currently being looped through		
		allocate(processDataHolder(size(particlePositions)/3, 3))
		processDataHolder = upperBound(1)*100
		
		if (dologging .and. rank == 0) print *, "--> Distributing particles among processors"
		
		!Loops through processors
		do counter2 = 1, nprocs
			!Stores how many particles are currently in this particles domain so it doesn't send too much data
			processTracker = 0	

			!This processors current lower and upper bound
			curr2(1:6) = processorRange(counter2, 1:6)

			!Loops through the 
			do counter1 = 1, size(particlePositions)/3 

				!When a particle is assigned a domain, it modifies the value so that it can be easily ignored
				if (particlePositions(counter1, 1) <= upperBound(1)) then

					!The current particle
					currParticle = particlePositions(counter1, 1:3)

					!If it's in the domain
					if (        currParticle(1) >= curr2(1) .and. currParticle(1) < curr2(4)&
						& .and. currParticle(2) >= curr2(2) .and. currParticle(2) < curr2(5)&
						& .and. currParticle(3) >= curr2(3) .and. currParticle(3) < curr2(6)) then
						!So it knows how many particles are stored in the current data holder, and adds it to said data holder
						processTracker = processTracker + 1
						processDataHolder(processTracker, 1:3) = currParticle(1:3)
							
						!Tagging the particle
						particlePositions(counter1, 1:3) = upperBound(1)*1000
					end if
				end if
			end do 

			!Now it knows the number of particles, its boundary, and its particle data, send it to the process

			!If it is itself, just assign the values
			if (counter2 == 1) then
				!How large it's original data set was
				arraySize = processTracker

				!Its upper and lower boundary
				lowerAndUpper(1:6) = curr2(1:6)

				allocate(sentProcessData(arraySize, 3))
				sentProcessData(1:arraySize, 1:3) = processDataHolder(1:arraySize, 1:3)
			else
				!numparticles
				call MPI_SEND(processTracker, 1, MPI_INTEGER, counter2-1, counter2-1, MPI_COMM_WORLD, ierr)

				!boundary
				call MPI_SEND(curr2(1:6), 6, MPI_DOUBLE, counter2-1, (counter2-1)*2, MPI_COMM_WORLD, ierr)

				!Particles
				call MPI_SEND(processDataHolder(1:processTracker, 1:3), 3*processTracker, MPI_DOUBLE, counter2-1, &
				&(counter2-1)*3, MPI_COMM_WORLD, ierr)
			end if
		end do

		deallocate(particlePositions)
		deallocate(processDataHolder)
		deallocate(processorRange)
	end if 

	!Now the other processes receive the corresponding data
	if (rank /= 0) then
		call MPI_RECV(arraySize, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, status1, ierr)
		call MPI_RECV(lowerAndUpper, 6, MPI_DOUBLE, 0, rank*2, MPI_COMM_WORLD, status1, ierr)

		allocate(sentProcessData(arraySize, 3))

		call MPI_RECV(sentProcessData, 3*arraySize, MPI_DOUBLE, 0, rank*3, MPI_COMM_WORLD, status1, ierr)
	end if


	!Send to all other processors the num of domains
	call MPI_BCAST(numDomains, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

	!And the cutoff incase it was updated
	call MPI_BCAST(cutoff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

	!The array which stores the processor ranks of the adjacent neighbours
	adjacents = getAdjacents(rank, numDomains)

	if (dologging .and. rank == 0) print *, "--> Distribution Complete"
	

	!Counts particles in it's own domain
	particlePairCount = 0

	if (dologging .and. rank == 0) print *, "--> Comparing particles in own process domain."	
	do counter2 = 1, arraySize
		do counter3 = counter2+1, arraySize
			 particlePairCount = particlePairCount + inRange(sentProcessData(counter2, 1:3),&
			 &sentProcessData(counter3, 1:3), cutoff, rank)
		end do
	end do

	call MPI_BARRIER(MPI_COMM_WORLD, ierr)

	if (dologging .and. rank == 0) print *, "--> Own particle comparisons finished."	

	!All data received will also store anything from neighbouring regions; initialise it with the sentProcessData
	allocate(allDataReceived(arraySize, 3))
	allDataReceived = sentProcessData

	!The current iteration
	counter2 = 0
 
	if (dologging .and. rank == 0) print *, "--> Starting iterations."	

	!The numbers represent the negative direction in the adjacents list
	do counter1 = 1, 5, 2

		counter2 = counter2 + 1 

		allocate(transferThese(size(allDataReceived)/3, 3))
		allocate(transferTheseAfter(size(allDataReceived)/3, 3)) 
		allocate(addAtEnd(size(allDataReceived)/3, 3))

		processTracker = 0 
		oppositeCount = 0 
		leftoverTracker = 0 

		!Seeing where the particles fall in the current axis in reference to sending to other processors
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

		!Sends/Receives the size of each relevant array to the neighbours
		numRecieved = 0
		CALL MPI_Sendrecv(processTracker, 1, MPI_INTEGER, adjacents(counter1), 2, &
           				     numRecieved, 1, MPI_INTEGER, adjacents(counter1+1), 2, &
           		 MPI_COMM_WORLD, status1, ierr)

		numRecievedAfter = 0

		CALL MPI_Sendrecv(oppositeCount, 1, MPI_INTEGER, adjacents(counter1+1), 3, &
           			   numRecievedAfter, 1, MPI_INTEGER, adjacents(counter1), 3, &
            		   MPI_COMM_WORLD, status1, ierr)

		!Allocates the arrays which will store the received particles temporarily
		allocate(receivedThese(numRecieved, 3))
		allocate(receivedTheseAfter(numRecievedAfter, 3))

		!Sending the particle data means first offsetting it from the relevant boundary vertex:
		!Eg: a particle in the very top corner can interact with the particle in the very bottom corner, but without modifying their positions they arent in cut off distance

		do counter3 = 1, processTracker
			transferThese(counter3, 1:3) = transferThese(counter3, 1:3) - lowerAndUpper(1:3)
		end do 

		if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		do counter3 = 1, oppositeCount
			transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) - currParticle
		end do 

		!Sends/Receives the particles to the processors
		CALL MPI_Sendrecv(transferThese(1:processTracker, 1:3), processTracker*3, MPI_DOUBLE, adjacents(counter1), 10, &
          					receivedThese(1:numRecieved, 1:3), numRecieved*3, MPI_DOUBLE, adjacents(counter1+1), 10, &
          			MPI_COMM_WORLD, status1, ierr)

		call MPI_Sendrecv(transferTheseAfter(1:oppositeCount, 1:3), oppositeCount*3, MPI_DOUBLE,  adjacents(counter1+1), 10, &
				    receivedTheseAfter(1:numRecievedAfter, 1:3), numRecievedAfter*3, MPI_DOUBLE,  adjacents(counter1), 10, &
				& MPI_COMM_WORLD, status1, ierr)

		!Adds the particle offsets received to the relevant boundary vertex to get the position of the particle in local coordinate space
		if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		do counter3 = 1, numRecieved
			receivedThese(counter3, 1:3) = receivedThese(counter3, 1:3) + currParticle
		end do

		do counter3 = 1, numRecievedAfter
			receivedTheseAfter(counter3, 1:3) = receivedTheseAfter(counter3, 1:3) + lowerAndUpper(1:3)
		end do

		!Undoing the offset calculation for its particles it already had so the positions aren't wrong for future iterations
		do counter3 = 1, processTracker
			transferThese(counter3, 1:3) = transferThese(counter3, 1:3) + lowerAndUpper(1:3)
		end do 

		if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3)]
		if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3)]
		if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6)]

		do counter3 = 1, oppositeCount
			transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) + currParticle
		end do 

		!Comparing the original particles it had to any newly received ones, but only in one direction to prevent double counting
		do i = 1, arraySize
			do j = 1, numRecieved
				 particlePairCount = particlePairCount + inRange(sentProcessData(i, 1:3),&
				 & receivedThese(j, 1:3), cutoff, rank)				
			end do
		end do


		!Transfering all particle data into allDataReceived, then deallocating any RAM-Style arrays to be used in the future iteration.
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

	deallocate(sentProcessData)
	deallocate(allDataReceived)

	!Send particle data back to the main process for final counting
	if (rank /= 0) then
		call MPI_SEND(particlePairCount, 1, MPI_INTEGER8, 0, 3, MPI_COMM_WORLD, ierr)
	else
		do i = 1, nprocs-1
			counter1 = 0
			call MPI_RECV(counter1, 1, MPI_INTEGER8, i, 3, MPI_COMM_WORLD, status1, ierr)
			particlePairCount = particlePairCount + counter1 
		end do
		print "(a46, f6.2, a3, i15, a5)", "<==> Total number of unique pairs with cutoff ", cutoff, " : ", particlePairCount, " <==>"
	end if


	call MPI_FINALIZE(ierr)

contains
	!Simply gets the modulus of the distance and compares against the cutoff, nothing fancy like PBC
	integer function inRange(p1, p2, cutoff, rank) result(count)
		double precision, dimension(1:3), intent(in) :: p1, p2
		double precision, dimension(1:3) :: temp
		double precision, intent(in):: cutoff
		integer, intent(in) :: rank

		temp = p2 - p1
		count = 0

		if (vectorMod(temp) < cutoff) count = 1
	end function inRange

	!For each process, it gets the adjacent processor ranks
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

		!This rather messy piece of code gets the upper and lower bound in each axis
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
	end function getAdjacents

	!Determines the lower and upper boundaries of each processor
	function getBounds(nprocs, lowerBound, upperBound, axisSplits) result(boundaryRanges)
		integer, intent(in) :: nprocs
		double precision, dimension(1:3), intent(in) :: lowerBound, upperBound
		double precision, dimension(nprocs,6) :: boundaryRanges
		double precision, dimension(1:3) :: difference, rangePerAxis
		integer, dimension(1:3), intent(in) :: axisSplits
		integer :: x, y, z, numIterations

		difference = upperBound - lowerBound
		
		rangePerAxis(1) = difference(1) / axisSplits(1)
		rangePerAxis(2) = difference(2) / axisSplits(2)
		rangePerAxis(3) = difference(3) / axisSplits(3)

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
	end function getBounds 

	!Determines the number of splits per axis to 'assign' to each process
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

		!breaks a given number into a product of prime factors
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

		!(Almost) Evenly distribute the factors among each major axis
		do i = 0, numFactors-1
			numDomains(1+mod(i, 3)) = numDomains(1+mod(i, 3)) * factors(i+1)
		end do 
	end function determineNoSplits

	!Function designed by george, making used of implicit save: used to check how long each section takes rather than the overall code when using time
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