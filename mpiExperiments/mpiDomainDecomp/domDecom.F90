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

	real, dimension(1:3) :: lowerBound, upperBound, curr1
	real, dimension(6) :: curr2
	integer :: counter1, counter2, counter3
	real :: cutoff

	real, dimension(1000,3) :: particlePositions

	integer :: particlePairCount
	integer :: totalSum, tempVar

	!mpi stuff
	integer :: ierr, rank, nprocs
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio, arraySize

	real, dimension(1:6) :: lowerAndUpper
	integer, dimension(1:6) :: adjacents
	real, dimension(:,:), allocatable :: processorRange, sentProcessData, processDataHolder !the 2nd : is going to be 6, lowerboundxyz, upperboundxyz
	real, dimension(:,:), allocatable :: lowerCoordsPerm, upperCoordsPerm, lowerCoordsChange, upperCoordsChange, leftovers !for use in each process to sort the data into the left/right, down/up, back/forward regions
	real, dimension(:,:), allocatable :: adParticles !for use with adjcanet process particle things
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
		processDataHolder = upperBound(1)*10 !arbitrary number so we know when its wrong
	
		!finding which particles go where and sending this to the various process'
		do counter2 = 1, nprocs !currentProcessr
			processTracker = 0
			processDataHolder = upperBound(1)*10
			curr2(1:6) = processorRange(counter2, 1:6)

			do counter1 = 1, size(particlePositions)/3 !currentParticle
				curr1 = particlePositions(counter1, 1:3)
				
				if (       curr1(1) >= curr2(1) .and. curr1(1) < curr2(4)&
					&.and. curr1(2) >= curr2(2) .and. curr1(2) < curr2(5)&
					&.and. curr1(3) >= curr2(3) .and. curr1(3) < curr2(6)) then
					processTracker = processTracker + 1
					processDataHolder(processTracker, 1:3) = curr1(1:3)
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

		!TODO: Allocate the correct amount of memory, be it through two messages, one with the size, one with data or just both

		! arraySize = processTracker(1) !
		! lowerAndUpper(1:6) =  processorRange(1, 1:6)

		! do counter1 = 1, nprocs-1 !represents the current process
		! 	call MPI_SEND(processTracker(counter1+1), 1, MPI_INTEGER, counter1, counter1, MPI_COMM_WORLD, ierr)
		! 	call MPI_SEND(processorRange(counter1+1, 1:6), 6, MPI_REAL, counter1, counter1*2, MPI_COMM_WORLD, ierr)
		! end do 
	end if 

	if (rank /= 0) then
		call MPI_RECV(arraySize, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, status1, ierr)
		call MPI_RECV(lowerAndUpper, 6, MPI_REAL, 0, rank*2, MPI_COMM_WORLD, status1, ierr)

		allocate(sentProcessData(arraySize, 3))

		call MPI_RECV(sentProcessData, 3*arraySize, MPI_REAL, 0, rank*3, MPI_COMM_WORLD, status1, ierr)
	end if

	!cant work with anything above 8 processors rn because one prioritses cutting the X dimension first whereas the other cuts the Z first

	adjacents = getAdjacantRegions(nprocs, rank)
	!print *, rank, ":", adjacents, "////", lowerAndUpper

	!first do itself
	!pair counting
	particlePairCount = 0
	do counter2 = 1, size(sentProcessData)/3
		do counter3 = counter2+1, size(sentProcessData)/3
			particlePairCount = particlePairCount + inRange(sentProcessData(counter2, 1:3),&
			&sentProcessData(counter3, 1:3), cutoff, rank)
		end do
	end do

	!print *, "Test 1:", rank, particlePairCount

	!print *, "Part 1", rank, particlePairCount
	!print *, rank, size(sentProcessData)/3, particlePairCount !sound as a cloud

	!now each process must send its particle data to the left and right
	!sending data to the left first, receiving from the right
	allocate(lowerCoordsPerm(size(sentProcessData)/3, 3))
	allocate(upperCoordsPerm(size(sentProcessData)/3, 3))
	allocate(leftovers(size(sentProcessData)/3, 3))
	
	lowerCoordsPerm = upperBound(1)*100
	upperCoordsPerm = upperBound(1)*100

	!this establishes where particles go in the new arrays
	curr2(1:3) = 0 !use curr2 to track the index of where to put coordinates in lower and upper coords arrays
	do counter1 = 1, size(sentProcessData)/3
		curr1(1:3) = sentProcessData(counter1, 1:3) !CURRENT PARTICLE, CHECK WHERE IT GOES
		!if (curr1(1) <= (lowerAndUpper(1) + ((lowerAndUpper(4)-lowerAndUpper(1))/2)) ) then !if the particles x coordinate is less than half way between the x boundaries
		if (curr1(1) <= (lowerAndUpper(1) + cutoff*(upperBound(1)-lowerBound(1)) ) ) then !if the particles x coordinate is less than half way between the x boundaries
			curr2(1) = curr2(1) + 1
			lowerCoordsPerm(int(curr2(1)), 1:3) = curr1(1:3)
		!else
		else if (curr1(1) <= (lowerAndUpper(4) - cutoff*(upperBound(1)-lowerBound(1)) ) ) then
			curr2(2) = curr2(2) + 1
			upperCoordsPerm(int(curr2(2)), 1:3) = curr1(1:3)
		else
			curr2(3) = curr2(3) + 1
			leftovers(int(curr2(3)), 1:3) = curr1(1:3)
		end if
	end do

	allocate(lowerCoordsChange(int(curr2(1)), 3))
	allocate(upperCoordsChange(int(curr2(2)), 3))

	!now you know how many each array are, 
	!create a copy into the changeable arrays

	!send to the left first the number of particles beign sent over
	call MPI_SEND(int(curr2(1)), 1, MPI_INTEGER, adjacents(1), rank, MPI_COMM_WORLD, ierr)

	!receive the size of the process to the right and allocate an array to store le data
	call MPI_RECV(counter1, 1, MPI_INTEGER, adjacents(2), adjacents(2), MPI_COMM_WORLD, status1, ierr)


	!now populate adParticles in a specific manner:
	!the particles that are being sent must be the offset from the lowest corner
	!when the particle data is received, in this case you add those vectors to this current boxes' highest x axis lowest everything else corner 

	do counter3 = 1, int(curr2(1))
		curr1(1:3) = lowerCoordsPerm(counter3, 1:3) !current particle
		lowerCoordsChange(counter3, 1:3) = curr1(1:3) - lowerAndUpper(1:3) 
	end do 
	
	!self explanatory
	do counter3 = 1, int(curr2(2))
		upperCoordsChange(counter3, 1:3) = upperCoordsPerm(counter3, 1:3)
	end do 

	allocate( adParticles((int(curr2(2)) + counter1), 3) )
	adParticles = upperBound(1)*100 !arbitrary

	adParticles(1:int(curr2(2)), 1:3) = upperCoordsChange(1:int(curr2(2)), 1:3)
	!adding the particles we have to the array where we will check coordinates

	! !send the data and receive the data
	call MPI_SEND(lowerCoordsChange(1:int(curr2(1)),1:3), int(curr2(1))*3, MPI_REAL, adjacents(1), rank,&
	& MPI_COMM_WORLD, ierr) !send the particles


	! !receive the data directly into adParticles
	call MPI_RECV(adParticles(int(curr2(2))+1:, 1:3), counter1*3, MPI_REAL, &
	& adjacents(2), adjacents(2), MPI_COMM_WORLD, status1, ierr)

	!now for those newly added things, add the offset
	do counter3 = int(curr2(2))+1, size(adParticles)/3
		adParticles(counter3 , 1:3) = adParticles(counter3 , 1:3) + [lowerAndUpper(4), lowerAndUpper(2), lowerAndUpper(3)]
	end do
	

	!///////////////////////////////
	!///////////////////////////////
	!///////////////////////////////	By this point, the adParticles array is fully populated with the particle data and now bounds must be calculated
	!///////////////////////////////	Between the first half of the particles adn the second half
	!///////////////////////////////

	!looping through the particles in a manner where we dont get duplicate results of self to self particles
	!no duplicate data on the first pass, though
	!if (rank > adjacents(2)) then !only the ones that have a higher rank shall do it in order to prevent duplicate data
		do counter2 = 1, int(curr2(2))
			do counter3 = int(curr2(2))+1, size(adParticles)/3
				particlePairCount = particlePairCount + inRange(adParticles(counter2, 1:3),&
				&adParticles(counter3, 1:3), cutoff, rank)
			end do
		end do		
	!end if		

	!print *, "Test 2:", rank, particlePairCount

	!You want to reuse adParticles for the next step, so restore and reallocate all this data in upperCoordsPerm as these new particles should remain eventually
	deallocate(upperCoordsPerm)
	allocate(upperCoordsPerm(size(adParticles)/3, 3))

	upperCoordsPerm = adParticles

	deallocate(adParticles)

	!also, reset lowerCoordinatesChange to match Perm
	lowerCoordsChange(1:int(curr2(1)), 1:3) = lowerCoordsPerm(1:int(curr2(1)), 1:3)


	!////////////////////////// Now do this all again but for left parickles

	!send the size of particles on the right
	!allocate the data for that
	!send that data
	!store it and add the offset
	

	!send to the right the number of particles beign sent over
	call MPI_SEND(int(curr2(2)), 1, MPI_INTEGER, adjacents(2), rank, MPI_COMM_WORLD, ierr)

	! !receive the size of the process to the left and allocate an array to store le data, stored in counter1
	call MPI_RECV(counter1, 1, MPI_INTEGER, adjacents(1), adjacents(1), MPI_COMM_WORLD, status1, ierr)	

	!the size of adParticles should now fit the data received from the left as whatever is on the left
	allocate(adParticles( int(curr2(1)) + counter1, 3))
	adParticles = upperBound(1)*100 !again, arbitrary

	!adding the known particles to the list
	adParticles(1:int(curr2(1)), 1:3) = lowerCoordsChange(1:int(curr2(1)), 1:3)

	!ensure you modify the upperCoordsChange data to be the vector not position
	do counter3 = 1, int(curr2(2))
				! if (rank == 1 .or. rank == 0) print *, rank, upperCoordsChange(counter3, 1:3), &
				! & upperCoordsChange(counter3, 1:3) - lowerAndUpper(1:3) 
		upperCoordsChange(counter3, 1) = upperCoordsChange(counter3, 1) - lowerAndUpper(4) !distance from the boundary which will connect to the neighbouring cell
		upperCoordsChange(counter3, 2) = upperCoordsChange(counter3, 2) - lowerAndUpper(2) 
		upperCoordsChange(counter3, 3) = upperCoordsChange(counter3, 3) - lowerAndUpper(3) 
	end do 
	
	!send the data and receive the data
	call MPI_SEND(upperCoordsChange(1:int(curr2(2)), 1:3), int(curr2(2))*3, MPI_REAL, adjacents(2), rank, &
	& MPI_COMM_WORLD, ierr) !send the particles


	! !receive the data directly into adParticles
	call MPI_RECV(adParticles(int(curr2(1))+1:, 1:3), counter1*3, MPI_REAL, &
	& adjacents(1), adjacents(1), MPI_COMM_WORLD, status1, ierr)

	!now for those newly added things, add the offset
	do counter3 = int(curr2(1))+1, size(adParticles)/3
		adParticles(counter3 , 1:3) = [ &
					& lowerAndUpper(1) + adParticles(counter3 , 1), &
					& lowerAndUpper(2) + adParticles(counter3 , 2), &
					& lowerAndUpper(3) + adParticles(counter3 , 3)&
				&]
	end do

	!now check pairs again
	if (rank <= adjacents(1)) then !if it received data from a box of a higher or equal rank, do something
		do counter2 = 1, int(curr2(1))
			do counter3 = int(curr2(1))+1, size(adParticles)/3
				particlePairCount = particlePairCount + inRange(adParticles(counter2, 1:3),&
				&adParticles(counter3, 1:3), cutoff, rank)
			end do
		end do		
	end if		
	!pairs added, now the new particles need to be stored in the perm variable and the rest can be freed, !!!!!!!!!!!!!!!!!!!!!!!!!X IS DONE

	!print *, "Test 3:", rank, particlePairCount

	deallocate(lowerCoordsPerm)
	allocate(lowerCoordsPerm(size(adParticles)/3, 3))

	lowerCoordsPerm = adParticles

	deallocate(adParticles)
	!now, upper and lower coords permanent should store all the particles of both adjacent regions

	!now, add these back into the main particle storing thign and repeat for alternate directions
	
	deallocate(sentProcessData)

	allocate(sentProcessData(size(lowerCoordsPerm)/3 + size(leftovers)/3 + size(upperCoordsPerm)/3, 3))

	sentProcessData(1:size(lowerCoordsPerm)/3, 1:3) = lowerCoordsPerm(1:size(lowerCoordsPerm)/3, 1:3)

	sentProcessData((size(lowerCoordsPerm)/3)+1:(size(lowerCoordsPerm)/3 + size(leftovers)/3), 1:3)&
	& = leftovers(1:size(upperCoordsPerm)/3, 1:3)

	sentProcessData((size(lowerCoordsPerm)/3 + size(leftovers)/3)+1:(size(lowerCoordsPerm)/3 + &
	& size(leftovers)/3 + size(upperCoordsPerm)/3), 1:3) = upperCoordsPerm(1:size(upperCoordsPerm)/3, 1:3)

	deallocate(lowerCoordsPerm)
	deallocate(lowerCoordsChange)
	deallocate(upperCoordsPerm)
	deallocate(upperCoordsChange)
	deallocate(leftovers)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! X Direction Done. Now, repeat the process in the Y direction.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(lowerCoordsPerm(size(sentProcessData)/3, 3))
	allocate(upperCoordsPerm(size(sentProcessData)/3, 3))
	allocate(leftovers(size(sentProcessData)/3, 3))

	lowerCoordsPerm = upperBound(1)*100
	upperCoordsPerm = upperBound(1)*100

	curr2(1:3) = 0 !use curr2 to track the index of where to put coordinates in lower and upper coords arrays

	! do counter1 = 1, size(sentProcessData)/3
	! 	curr1(1:3) = sentProcessData(counter1, 1:3) !CURRENT PARTICLE, CHECK WHERE IT GOES
	! 	if (curr1(2) <= (lowerAndUpper(2) + ((lowerAndUpper(5)-lowerAndUpper(2))/2)) ) then !if the particles x coordinate is less than half way between the x boundaries
	! 		curr2(1) = curr2(1) + 1
	! 		lowerCoordsPerm(int(curr2(1)), 1:3) = curr1(1:3)
	! 	else
	! 		curr2(2) = curr2(2) + 1
	! 		upperCoordsPerm(int(curr2(2)), 1:3) = curr1(1:3)
	! 	end if
	! end do
	
	do counter1 = 1, size(sentProcessData)/3
		curr1(1:3) = sentProcessData(counter1, 1:3) !CURRENT PARTICLE, CHECK WHERE IT GOES
		!if (curr1(1) <= (lowerAndUpper(1) + ((lowerAndUpper(4)-lowerAndUpper(1))/2)) ) then !if the particles x coordinate is less than half way between the x boundaries
		if (curr1(2) <= (lowerAndUpper(2) + cutoff*(upperBound(2)-lowerBound(2)) ) ) then !if the particles x coordinate is less than half way between the x boundaries
			curr2(1) = curr2(1) + 1
			lowerCoordsPerm(int(curr2(1)), 1:3) = curr1(1:3)
		!else
		else if (curr1(2) <= (lowerAndUpper(5) - cutoff*(upperBound(2)-lowerBound(2)) ) ) then
			curr2(2) = curr2(2) + 1
			upperCoordsPerm(int(curr2(2)), 1:3) = curr1(1:3)
		else
			curr2(3) = curr2(3) + 1
			leftovers(int(curr2(3)), 1:3) = curr1(1:3)
		end if
	end do

	allocate(lowerCoordsChange(int(curr2(1)), 3))
	allocate(upperCoordsChange(int(curr2(2)), 3))

	!Sending beneath first
	call MPI_SEND(int(curr2(1)), 1, MPI_INTEGER, adjacents(3), rank, MPI_COMM_WORLD, ierr)

	!receive the size of the process to the right and allocate an array to store le data
	call MPI_RECV(counter1, 1, MPI_INTEGER, adjacents(4), adjacents(4), MPI_COMM_WORLD, status1, ierr)

	!this is finding the vector. You are adding the distance to the top of each box essentuially, so you just find this distance
	do counter3 = 1, int(curr2(1))
		curr1(1:3) = lowerCoordsPerm(counter3, 1:3) !current particle
		lowerCoordsChange(counter3, 1:3) = curr1(1:3) - lowerAndUpper(1:3) 
	end do 
	
	!self explanatory
	do counter3 = 1, int(curr2(2))
		upperCoordsChange(counter3, 1:3) = upperCoordsPerm(counter3, 1:3)
	end do 

	!so it can store all data coming in and here already
	allocate( adParticles((int(curr2(2)) + counter1), 3) )
	adParticles = upperBound(1)*100 !arbitrary

	adParticles(1:int(curr2(2)), 1:3) = upperCoordsChange(1:int(curr2(2)), 1:3)
	!adding the particles we have to the array where we will check coordinates

	! !send the data and receive the data
	call MPI_SEND(lowerCoordsChange(1:int(curr2(1)),1:3), int(curr2(1))*3, MPI_REAL, adjacents(3), rank,&
	& MPI_COMM_WORLD, ierr) !send the particles


	! !receive the data directly into adParticles
	call MPI_RECV(adParticles(int(curr2(2))+1:, 1:3), counter1*3, MPI_REAL, &
	& adjacents(4), adjacents(4), MPI_COMM_WORLD, status1, ierr)

	!now for those newly added things, add the offset
	do counter3 = int(curr2(2))+1, size(adParticles)/3
		adParticles(counter3 , 1:3) = adParticles(counter3 , 1:3) + [lowerAndUpper(1), lowerAndUpper(5), lowerAndUpper(3)] !the 2nd y this time not x
	end do

	!now that we already have some duplicate data we need to check if the particle you are currently dealing with lies within your boundaries. If it does, check, otherwise do not
	do counter2 = 1, int(curr2(2))
		do counter3 = int(curr2(2))+1, size(adParticles)/3
			if (adParticles(counter2, 1) >= lowerAndUpper(1) .and. adParticles(counter2, 1) <= lowerAndUpper(4)) then !logic to prevent duplicate data, revisit if needed
				particlePairCount = particlePairCount + inRange(adParticles(counter2, 1:3),&
				&adParticles(counter3, 1:3), cutoff, rank)
			end if
		end do
	end do		


	!You want to reuse adParticles for the next step, so restore and reallocate all this data in upperCoordsPerm as these new particles should remain eventually
	deallocate(upperCoordsPerm)
	allocate(upperCoordsPerm(size(adParticles)/3, 3))

	upperCoordsPerm = adParticles

	deallocate(adParticles)

	!also, reset lowerCoordinatesChange to match Perm
	lowerCoordsChange(1:int(curr2(1)), 1:3) = lowerCoordsPerm(1:int(curr2(1)), 1:3)


	!////////////////////////// Now do this all again but for left parickles

	!send the size of particles on the right
	!allocate the data for that
	!send that data
	!store it and add the offset
	

	!send to the right the number of particles beign sent over
	call MPI_SEND(int(curr2(2)), 1, MPI_INTEGER, adjacents(4), rank, MPI_COMM_WORLD, ierr)

	! !receive the size of the process to the left and allocate an array to store le data, stored in counter1
	call MPI_RECV(counter1, 1, MPI_INTEGER, adjacents(3), adjacents(3), MPI_COMM_WORLD, status1, ierr)	

	!the size of adParticles should now fit the data received from the left as whatever is on the left
	allocate(adParticles( int(curr2(1)) + counter1, 3))
	adParticles = upperBound(1)*100 !again, arbitrary

	!adding the known particles to the list
	adParticles(1:int(curr2(1)), 1:3) = lowerCoordsChange(1:int(curr2(1)), 1:3)

	!ensure you modify the upperCoordsChange data to be the vector not position
	do counter3 = 1, int(curr2(2))
		upperCoordsChange(counter3, 1) = upperCoordsChange(counter3, 1) - lowerAndUpper(1) !distance from the boundary which will connect to the neighbouring cell
		upperCoordsChange(counter3, 2) = upperCoordsChange(counter3, 2) - lowerAndUpper(5) 
		upperCoordsChange(counter3, 3) = upperCoordsChange(counter3, 3) - lowerAndUpper(3) 
	end do 
	
	!send the data and receive the data
	call MPI_SEND(upperCoordsChange(1:int(curr2(2)), 1:3), int(curr2(2))*3, MPI_REAL, adjacents(4), rank, &
	& MPI_COMM_WORLD, ierr) !send the particles

	! !receive the data directly into adParticles
	call MPI_RECV(adParticles(int(curr2(1))+1:, 1:3), counter1*3, MPI_REAL, &
	& adjacents(3), adjacents(3), MPI_COMM_WORLD, status1, ierr)

	!now for those newly added things, add the offset
	do counter3 = int(curr2(1))+1, size(adParticles)/3
		adParticles(counter3 , 1:3) = [ &
					& lowerAndUpper(1) + adParticles(counter3 , 1), &
					& lowerAndUpper(2) + adParticles(counter3 , 2), &
					& lowerAndUpper(3) + adParticles(counter3 , 3)&
				&]
	end do

	!now check pairs again
	if (rank <= adjacents(3)) then !if it received data from a box of a higher or equal rank, do something
		do counter2 = 1, int(curr2(1))
			do counter3 = int(curr2(1))+1, size(adParticles)/3
				if (adParticles(counter2, 1) >= lowerAndUpper(1) .and. adParticles(counter2, 1) <= lowerAndUpper(4)) then !logic to prevent duplicate data, revisit if needed
					particlePairCount = particlePairCount + inRange(adParticles(counter2, 1:3),&
					&adParticles(counter3, 1:3), cutoff, rank)
				end if
			end do
		end do		
	end if		
	!pairs added, now the new particles need to be stored in the perm variable and the rest can be freed, !!!!!!!!!!!!!!!!!!!!!!!!!X IS DONE

	!print *, "Test 3:", rank, particlePairCount

	deallocate(lowerCoordsPerm)
	allocate(lowerCoordsPerm(size(adParticles)/3, 3))

	lowerCoordsPerm = adParticles

	deallocate(adParticles)
	!now, upper and lower coords permanent should store all the particles of both adjacent regions

	!now, add these back into the main particle storing thign and repeat for alternate directions
	
	deallocate(sentProcessData)

	allocate(sentProcessData(size(lowerCoordsPerm)/3 + size(leftovers)/3 + size(upperCoordsPerm)/3, 3))

	sentProcessData(1:size(lowerCoordsPerm)/3, 1:3) = lowerCoordsPerm(1:size(lowerCoordsPerm)/3, 1:3)

	sentProcessData((size(lowerCoordsPerm)/3)+1:(size(lowerCoordsPerm)/3 + size(leftovers)/3), 1:3)&
	& = leftovers(1:size(upperCoordsPerm)/3, 1:3)

	sentProcessData((size(lowerCoordsPerm)/3 + size(leftovers)/3)+1:(size(lowerCoordsPerm)/3 + &
	& size(leftovers)/3 + size(upperCoordsPerm)/3), 1:3) = upperCoordsPerm(1:size(upperCoordsPerm)/3, 1:3)

	deallocate(lowerCoordsPerm)
	deallocate(lowerCoordsChange)
	deallocate(upperCoordsPerm)
	deallocate(upperCoordsChange)
	deallocate(leftovers)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !now for the Z
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(lowerCoordsPerm(size(sentProcessData)/3, 3))
	allocate(upperCoordsPerm(size(sentProcessData)/3, 3))
	allocate(leftovers(size(sentProcessData)/3, 3))

	lowerCoordsPerm = upperBound(1)*100
	upperCoordsPerm = upperBound(1)*100

	curr2(1:3) = 0 !use curr2 to track the index of where to put coordinates in lower and upper coords arrays

	! do counter1 = 1, size(sentProcessData)/3
	! 	curr1(1:3) = sentProcessData(counter1, 1:3) !CURRENT PARTICLE, CHECK WHERE IT GOES
	! 	if (curr1(3) <= (lowerAndUpper(3) + ((lowerAndUpper(6)-lowerAndUpper(3))/2)) ) then !if the particles x coordinate is less than half way between the x boundaries
	! 		curr2(1) = curr2(1) + 1
	! 		lowerCoordsPerm(int(curr2(1)), 1:3) = curr1(1:3)
	! 	else
	! 		curr2(2) = curr2(2) + 1
	! 		upperCoordsPerm(int(curr2(2)), 1:3) = curr1(1:3)
	! 	end if
	! end do

	do counter1 = 1, size(sentProcessData)/3
		curr1(1:3) = sentProcessData(counter1, 1:3) !CURRENT PARTICLE, CHECK WHERE IT GOES
		!if (curr1(1) <= (lowerAndUpper(1) + ((lowerAndUpper(4)-lowerAndUpper(1))/2)) ) then !if the particles x coordinate is less than half way between the x boundaries
		if (curr1(3) <= (lowerAndUpper(3) + cutoff*(upperBound(3)-lowerBound(3)) ) ) then !if the particles x coordinate is less than half way between the x boundaries
			curr2(1) = curr2(1) + 1
			lowerCoordsPerm(int(curr2(1)), 1:3) = curr1(1:3)
		!else
		else if (curr1(3) <= (lowerAndUpper(6) - cutoff*(upperBound(3)-lowerBound(3)) ) ) then
			curr2(2) = curr2(2) + 1
			upperCoordsPerm(int(curr2(2)), 1:3) = curr1(1:3)
		else
			curr2(3) = curr2(3) + 1
			leftovers(int(curr2(3)), 1:3) = curr1(1:3)
		end if
	end do

	allocate(lowerCoordsChange(int(curr2(1)), 3))
	allocate(upperCoordsChange(int(curr2(2)), 3))

	!Sending beneath first
	call MPI_SEND(int(curr2(1)), 1, MPI_INTEGER, adjacents(5), rank, MPI_COMM_WORLD, ierr)

	!receive the size of the process to the right and allocate an array to store le data
	call MPI_RECV(counter1, 1, MPI_INTEGER, adjacents(6), adjacents(6), MPI_COMM_WORLD, status1, ierr)

	!this is finding the vector. You are adding the distance to the top of each box essentuially, so you just find this distance
	do counter3 = 1, int(curr2(1))
		curr1(1:3) = lowerCoordsPerm(counter3, 1:3) !current particle
		lowerCoordsChange(counter3, 1:3) = curr1(1:3) - lowerAndUpper(1:3) 
	end do 
	
	!self explanatory
	do counter3 = 1, int(curr2(2))
		upperCoordsChange(counter3, 1:3) = upperCoordsPerm(counter3, 1:3)
	end do 

	!so it can store all data coming in and here already
	allocate( adParticles((int(curr2(2)) + counter1), 3) )
	adParticles = upperBound(1)*100 !arbitrary

	adParticles(1:int(curr2(2)), 1:3) = upperCoordsChange(1:int(curr2(2)), 1:3)
	!adding the particles we have to the array where we will check coordinates

	! !send the data and receive the data
	call MPI_SEND(lowerCoordsChange(1:int(curr2(1)),1:3), int(curr2(1))*3, MPI_REAL, adjacents(5), rank,&
	& MPI_COMM_WORLD, ierr) !send the particles


	! !receive the data directly into adParticles
	call MPI_RECV(adParticles(int(curr2(2))+1:, 1:3), counter1*3, MPI_REAL, &
	& adjacents(6), adjacents(6), MPI_COMM_WORLD, status1, ierr)

	!now for those newly added things, add the offset
	do counter3 = int(curr2(2))+1, size(adParticles)/3
		adParticles(counter3 , 1:3) = adParticles(counter3 , 1:3) + [lowerAndUpper(1), lowerAndUpper(2), lowerAndUpper(5)] !the 2nd z
	end do

	!now that we already have some duplicate data we need to check if the particle you are currently dealing with lies within your boundaries. If it does, check, otherwise do not
	do counter2 = 1, int(curr2(2))
		do counter3 = int(curr2(2))+1, size(adParticles)/3
			if (adParticles(counter2, 1) >= lowerAndUpper(1) .and. adParticles(counter2, 1) <= lowerAndUpper(4)&
			   & .and. adParticles(counter2, 2) >= lowerAndUpper(2) .and. adParticles(counter2, 2) <= lowerAndUpper(5)) then !logic to prevent duplicate data, revisit if needed
				particlePairCount = particlePairCount + inRange(adParticles(counter2, 1:3),&
				&adParticles(counter3, 1:3), cutoff, rank)
			end if
		end do
	end do		


	!You want to reuse adParticles for the next step, so restore and reallocate all this data in upperCoordsPerm as these new particles should remain eventually
	deallocate(upperCoordsPerm)
	allocate(upperCoordsPerm(size(adParticles)/3, 3))

	upperCoordsPerm = adParticles

	deallocate(adParticles)

	!also, reset lowerCoordinatesChange to match Perm
	lowerCoordsChange(1:int(curr2(1)), 1:3) = lowerCoordsPerm(1:int(curr2(1)), 1:3)

	call MPI_SEND(int(curr2(2)), 1, MPI_INTEGER, adjacents(6), rank, MPI_COMM_WORLD, ierr)

	! !receive the size of the process to the left and allocate an array to store le data, stored in counter1
	call MPI_RECV(counter1, 1, MPI_INTEGER, adjacents(5), adjacents(5), MPI_COMM_WORLD, status1, ierr)	

	!the size of adParticles should now fit the data received from the left as whatever is on the left
	allocate(adParticles( int(curr2(1)) + counter1, 3))
	adParticles = upperBound(1)*100 !again, arbitrary

	!adding the known particles to the list
	adParticles(1:int(curr2(1)), 1:3) = lowerCoordsChange(1:int(curr2(1)), 1:3)

	!ensure you modify the upperCoordsChange data to be the vector not position
	do counter3 = 1, int(curr2(2))
		upperCoordsChange(counter3, 1) = upperCoordsChange(counter3, 1) - lowerAndUpper(1) 
		upperCoordsChange(counter3, 2) = upperCoordsChange(counter3, 2) - lowerAndUpper(2) 
		upperCoordsChange(counter3, 3) = upperCoordsChange(counter3, 3) - lowerAndUpper(6) 
	end do 
	
	!send the data and receive the data
	call MPI_SEND(upperCoordsChange(1:int(curr2(2)), 1:3), int(curr2(2))*3, MPI_REAL, adjacents(6), rank, &
	& MPI_COMM_WORLD, ierr) !send the particles

	! !receive the data directly into adParticles
	call MPI_RECV(adParticles(int(curr2(1))+1:, 1:3), counter1*3, MPI_REAL, &
	& adjacents(5), adjacents(5), MPI_COMM_WORLD, status1, ierr)

	!now for those newly added things, add the offset
	do counter3 = int(curr2(1))+1, size(adParticles)/3
		adParticles(counter3 , 1:3) = [ &
					& lowerAndUpper(1) + adParticles(counter3 , 1), &
					& lowerAndUpper(2) + adParticles(counter3 , 2), &
					& lowerAndUpper(3) + adParticles(counter3 , 3)&
				&]
	end do

	!now check pairs again
	if (rank <= adjacents(5)) then !if it received data from a box of a higher or equal rank, do something
		do counter2 = 1, int(curr2(1))
			do counter3 = int(curr2(1))+1, size(adParticles)/3
				if (adParticles(counter2, 1) >= lowerAndUpper(1) .and. adParticles(counter2, 1) <= lowerAndUpper(4)&
				& .and. adParticles(counter2, 2) >= lowerAndUpper(2) .and. adParticles(counter2, 2) <= lowerAndUpper(5)) then !logic to prevent duplicate data, revisit if needed
					particlePairCount = particlePairCount + inRange(adParticles(counter2, 1:3),&
					&adParticles(counter3, 1:3), cutoff, rank)
				end if
			end do
		end do		
	end if		
	!pairs added, now the new particles need to be stored in the perm variable and the rest can be freed, !!!!!!!!!!!!!!!!!!!!!!!!!X IS DONE

	!print *, "Test 3:", rank, particlePairCount

	deallocate(lowerCoordsPerm)
	allocate(lowerCoordsPerm(size(adParticles)/3, 3))

	lowerCoordsPerm = adParticles

	deallocate(adParticles)
	!now, upper and lower coords permanent should store all the particles of both adjacent regions

	!now, add these back into the main particle storing thign and repeat for alternate directions
	
	deallocate(sentProcessData)


	allocate(sentProcessData(size(lowerCoordsPerm)/3 + size(leftovers)/3 + size(upperCoordsPerm)/3, 3))

	sentProcessData(1:size(lowerCoordsPerm)/3, 1:3) = lowerCoordsPerm(1:size(lowerCoordsPerm)/3, 1:3)

	sentProcessData((size(lowerCoordsPerm)/3)+1:(size(lowerCoordsPerm)/3 + size(leftovers)/3), 1:3)&
	& = leftovers(1:size(upperCoordsPerm)/3, 1:3)

	sentProcessData((size(lowerCoordsPerm)/3 + size(leftovers)/3)+1:(size(lowerCoordsPerm)/3 + &
	& size(leftovers)/3 + size(upperCoordsPerm)/3), 1:3) = upperCoordsPerm(1:size(upperCoordsPerm)/3, 1:3)

	deallocate(lowerCoordsPerm)
	deallocate(lowerCoordsChange)
	deallocate(upperCoordsPerm)
	deallocate(upperCoordsChange)
	deallocate(leftovers)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if ( rank /= 0) then
		call MPI_SEND(particlePairCount, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, ierr)
	else
		do counter1 = 1, nprocs-1
			counter2 = 0
			call MPI_RECV(counter2, 1, MPI_INTEGER, counter1, counter1, MPI_COMM_WORLD, status1, ierr)	
			particlePairCount = particlePairCount + counter2
		end do
		
		print *, "Total number of pairs:", particlePairCount
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