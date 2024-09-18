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

	real, dimension(10000,3) :: particlePositions

	integer, dimension(10000) :: particlePairCount
	integer :: totalSum, tempVar

	!mpi stuff
	integer :: ierr, rank, nprocs
    integer, dimension(MPI_STATUS_SIZE) :: status1
	integer :: particleProcRatio

	integer, dimension(:,:,:), allocatable :: processDataHolder, sentProcessData !the 3rd : is going to be 3 regardless in both
	real, dimension(:,:), allocatable :: processorRange !the 2nd : is going to be 6, lowerboundxyz, upperboundxyz
	integer, dimension(:) , allocatable :: processTracker
	lowerBound = [0.0, 0.0, 0.0]  !
	upperBound = [1.0, 1.0, 1.0] !
	cutoff = 0.49
	totalSum = 0
	tempVar = 0

	do counter1 = 1, size(particlePositions)/3
		particlePairCount(counter1) = 0
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Receiving position data from dat file
	!!! And populating position/velocity arrays 

	!Positions
	!Not strictly necessary
	!9 format(f20.17, 4x, f20.17, 4x, f20.17)

	!open(11, file="pData.dat", status="old")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!initlaise mpi 
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
		allocate(processDataHolder(nprocs,size(particlePositions)/3, 3))
		allocate(processTracker(nprocs))
		processDataHolder = 0
		processTracker = 0
	
		do counter1 = 1, size(particlePositions)/3
			curr1 = particlePositions(counter1, 1:3)
			do counter2 = 1, nprocs
				curr2(1:6) = processorRange(counter2, 1:6)
				if (       curr1(1) >= curr2(1) .and. curr1(1) < curr2(4)&
					&.and. curr1(2) >= curr2(2) .and. curr1(2) < curr2(5)&
					&.and. curr1(3) >= curr2(3) .and. curr1(3) < curr2(6)) then
					processTracker(counter2) = processTracker(counter2) + 1
					processDataHolder(counter2, processTracker(counter2), 1:3) = curr1
				end if
			end do 
		end do

		do counter1 = 1, nprocs
			print*, counter1, processTracker(counter1), processorRange(counter1, 1:6)
		end do 

		!TODO: Send (or broadcast) messages to the needed processors
		!TODO: Allocate the correct amount of memory, be it through two messages, one with the size, one with data or just both
		deallocate(processorRange)
		!see what the dimensions
	end if 

	call MPI_FINALIZE(ierr)

contains
	!its important to note that each process
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