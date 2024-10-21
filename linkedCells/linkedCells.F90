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
    character :: logging, writeToFiles, logTime
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

    print *, "-> Write Particles to a File (Y/N, Case Sensitive): "
    read(*,*) writeToFiles

    print *, "-> Log Time Results (Y/N, Case Sensitive): "
    read(*,*) logTime

    ! The main program uses an upper and lower bound, however this wants a minimum and a dimension as it makes changing the regions better
    axisDimension = axisDimension - minBoundary
    seed = c1       ! Not directly stored into seed as it will need 8 inputs every time which is a headache and unecessary for this program

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
    !SEED
    write(11, *) seed(1)
    !WritenInAfile
    write(11, *) writeToFiles
    !log time
    write(11, *) logTime

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
    
    ret_mod = v1(1)**2 + v1(2)**2 + v1(3)**2
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
  double precision, dimension(6) :: curr2     !RAM-functionality variables
  double precision, dimension(4) :: currParticle   !RAM-functionality variables
  integer, dimension(1:3) :: numDomains
  double precision, dimension(1:6) :: lowerAndUpper
  integer, dimension(1:6) :: adjacents            !What processors are Left,Right,Down,Up,Forward,Backward

  !Counter variables
  integer :: counter1, counter2, counter3, i, j

  !Output-Relevant variables
  double precision :: cutoff
  integer(kind=int64) :: particlePairCount, rootCount
  integer :: processTracker
  logical :: success, dologging
  character :: generateNewFiles, logging, writeToFiles, logTime
  integer :: particleProcRatio, arraySize, oppositeCount, leftoverTracker

  !Output-Relevant, variable sized arrays
  double precision, dimension(:,:), allocatable :: particlePositions                    !all particle positions
  double precision, dimension(:,:), allocatable :: processorRange, sentProcessData, processDataHolder
  double precision, dimension(:,:), allocatable :: transferThese, receivedThese
  double precision, dimension(:,:), allocatable :: transferTheseAfter, receivedTheseAfter, allDataReceived 
  double precision, dimension(:,:), allocatable :: addAtEnd
  real(kind=8), dimension(18) :: times, rootTimes
  integer, dimension(8) :: seed

  !mpi stuff
  integer :: ierr, rank, nprocs, numParticles, numRecieved, numRecievedAfter
  integer, dimension(MPI_STATUS_SIZE) :: status1

  logical :: testing = .true. !If true, it skips the USER input section 

  !Init for MPI
  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)  

  if (rank == 0) print *, "[Time] Init:", elapsed_time(), timestamp()

  ! Confirming if new data should be created or not. Checks if relevant files exist, if not program ends with message
  if (rank == 0 .and. .not.  testing) then
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

  if (rank == 0) i = timestamp()
  times = 0.00000000000000d0

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  open(10, file="config.txt", status="old")
  
  read(10, *) numParticles
  read(10, *) cutoff
  read(10, *) lowerBound
  read(10, *) upperBound
  read(10, *) logging
  read(10, *) i !seed
  read(10, *) writeToFiles
  read(10, *) logTime

  close(10)
  ! Master-process only things
  if (rank == 0) then
    seed = i
    call random_seed(put=seed)
    !Start the main program Timer

    !Initialiser for logging
    dologging = .false.
    if (logging == "Y") dologging = .true.

    !Priting config parameters
    print *, ""
    print *, "[OUTP] Program Config Parameters [OUTP]"
    print *, " - Cutoff:", cutoff
    print *, " - Number of Particles:", numParticles
    print *, " - Lower Bound:", lowerBound
    print *, " - Upper Bound:", upperBound
    print *, " - Logging Enabled: ", dologging
    print *, " - Saved in File: ", writeToFiles
    print *, " - Logging Time: ", logTime
    print *, "[OUTP] End of List [OUTP]"
    print "(//a28/)", " [OUTP] Program Start [OUTP]"
    if (logTime == "Y") print *, "[TIME] Getting times..."

    if (dologging .and. rank == 0) print *, "[OUTP] Reading file data"
    !Want to be able to gloss over the time it takes to do anything relating to changing the particles or config file
    i = elapsed_time()

    allocate(particlePositions(numParticles, 4))
    
    ! All particle positions, only known for this process at this time
    inquire(file="coordinates.txt", exist=success)
    if (success .and. writeToFiles == "Y") then
      open(11, file="coordinates.txt", status="old")

      do counter1 = 1, numParticles
        read(11,*) particlePositions(counter1, 1), particlePositions(counter1, 2), particlePositions(counter1, 3)
        particlePositions(counter1, 4) = counter1
      end do 

      close(11)
    else 
      call random_number(particlePositions)
      !Formatting the particles to be distributed evenly across the given upper and lower bound
      do counter1 = 1, numParticles
        particlePositions(counter1, 1:3) = [&
          & lowerBound(1) + particlePositions(counter1, 1)*(upperBound(1) - lowerBound(1)), &
          & lowerBound(2) + particlePositions(counter1, 2)*(upperBound(2) - lowerBound(2)), &
          & lowerBound(3) + particlePositions(counter1, 3)*(upperBound(3) - lowerBound(3))  &
        &]
      end do
      !Generating the particles is much faster
    end if 

    if (dologging .and. rank == 0) print *, "[OUTP] Particle Data Generated"
    
    !if (rank == 0 .and. logTime == "Y") print *, "[TIME] Particle Generation:", elapsed_time()
    times(1) = elapsed_time()

    ! Number of slices in each axis to distribute processors
    numDomains = determineNoSplits(nprocs)

    ! In case an improper cutoff is put in, checks it
    cutoff = min( min(cutoff, (upperBound(1)-lowerBound(1))/ (numDomains(1)*2)), &
      min( (upperBound(2)-lowerBound(2))/(numDomains(2)*2), (upperBound(3)-lowerBound(3))/(numDomains(3)*2)) )

    12 format(" [OUTP] Space divided across:", i4, " processors, Number of domains per axis:", i4, ",", i4, "," i4)
    if (dologging .and. rank == 0) print 12, nprocs, numDomains(1), numDomains(2), numDomains(3)

    !if (rank == 0 .and. logTime == "Y") print *, "[Time] Divide Axis:", elapsed_time()
    times(2) = elapsed_time()
    !The upper and lower bound of each process
    allocate(processorRange(nprocs,6))
    processorRange = getBounds(nprocs, lowerBound, upperBound, numDomains)
    
    if (dologging .and. rank == 0) print *, "[OUTP] Distributing particles among processors"
    i = elapsed_time()

    do counter2 = 1, nprocs
      curr2(1:6) = processorRange(counter2, 1:6)
      if (counter2 == 1) then
        !Its upper and lower boundary
        lowerAndUpper(1:6) = curr2(1:6)
      else
        !boundary
        call MPI_SEND(curr2(1:6), 6, MPI_DOUBLE, counter2-1, (counter2-1)*2, MPI_COMM_WORLD, ierr)
      end if
    end do

  end if 

  !!TRYING TO MAKE ALL PROCESSORS READ PARTICLES AT THE SAME TIME 
  !Send to all other processors the num of domains
  call MPI_BCAST(numParticles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  if (rank /= 0) allocate(particlePositions(numParticles, 4))
  call MPI_BCAST(particlePositions, numParticles*4, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

  if (rank /= 0) then
    call MPI_RECV(lowerAndUpper, 6, MPI_DOUBLE, 0, rank*2, MPI_COMM_WORLD, status1, ierr)
  end if
    
  allocate(processDataHolder(size(particlePositions)/4, 4))
  processDataHolder = upperBound(1)*100

  processTracker = 0  

  do counter1 = 1, size(particlePositions)/4 
    currParticle = particlePositions(counter1, 1:4)

    if (      currParticle(1) >= lowerAndUpper(1) .and. currParticle(1) < lowerAndUpper(4)&
      & .and. currParticle(2) >= lowerAndUpper(2) .and. currParticle(2) < lowerAndUpper(5)&
      & .and. currParticle(3) >= lowerAndUpper(3) .and. currParticle(3) < lowerAndUpper(6)) then
      processTracker = processTracker + 1
      processDataHolder(processTracker, 1:4) = currParticle(1:4)
    end if
  end do 

  allocate(sentProcessData(processTracker, 4))
  arraySize = processTracker

  sentProcessData(1:arraySize, 1:4) = processDataHolder(1:arraySize, 1:4)
  deallocate(processDataHolder)

  !Send to all other processors the num of domains
  call MPI_BCAST(numDomains, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !And the cutoff incase it was updated
  call MPI_BCAST(cutoff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

  if (dologging .and. rank == 0) print *, "[OUTP] Distribution Complete"

  times(3) = elapsed_time()

  adjacents = getAdjacents(rank, numDomains)
  
  times(4) = elapsed_time()

  particlePairCount = 0

  if (dologging .and. rank == 0) print *, "[OUTP] Comparing particles in own process domain."  
  do counter2 = 1, arraySize
    do counter3 = counter2+1, arraySize
       particlePairCount = particlePairCount + inRange(sentProcessData(counter2, 1:3),&
       &sentProcessData(counter3, 1:3), cutoff, rank)
    end do
  end do

  !if (rank == 0 .and. logTime == "Y") print *, "[Time] Comparing Own Particles:", elapsed_time()
  times(5) = elapsed_time()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if (dologging .and. rank == 0) print *, "[OUTP] Own particle comparisons finished."  

  !All data received will also store anything from neighbouring regions; initialise it with the sentProcessData
  allocate(allDataReceived(arraySize, 4))
  allDataReceived = sentProcessData

  !The current iteration
  counter2 = 0
 
  if (dologging .and. rank == 0) print *, "[OUTP] Starting iterations."  
  !The numbers represent the negative direction in the adjacents list
  do counter1 = 1, 5, 2
    
    counter2 = counter2 + 1 
    i = elapsed_time() !so it resets
    allocate(transferThese(size(allDataReceived)/4, 4))
    allocate(transferTheseAfter(size(allDataReceived)/4, 4)) 
    allocate(addAtEnd(size(allDataReceived)/4, 4))

    processTracker = 0 
    oppositeCount = 0 
    leftoverTracker = 0 

    !Seeing where the particles fall in the current axis in reference to sending to other processors
    do counter3 = 1, size(allDataReceived)/4
      currParticle = allDataReceived(counter3, 1:4)
      if ( currParticle(counter2) < lowerAndUpper(counter2)+cutoff ) then
        processTracker = processTracker + 1
        transferThese(processTracker, 1:4) = currParticle(1:4)
      else if ( currParticle(counter2) >= lowerAndUpper(counter2+3)-cutoff) then
        oppositeCount = oppositeCount + 1
        transferTheseAfter(oppositeCount, 1:4) = currParticle(1:4)
      else
        leftoverTracker = leftoverTracker + 1
        addAtEnd(leftoverTracker, 1:4) = currParticle(1:4)
      end if
    end do

    !if (rank == 0 .and. logTime == "Y") print "(a18,i1,a30,f20.16)", " [TIME] Iteration ", counter2, &
    !& " Particle Region Calculation: ", elapsed_time()
    times(5 + 4*(counter2-1) + 1) = elapsed_time()
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
    allocate(receivedThese(numRecieved, 4))
    allocate(receivedTheseAfter(numRecievedAfter, 4))

    !Sending the particle data means first offsetting it from the relevant boundary vertex:
    !Eg: a particle in the very top corner can interact with the particle in the very bottom corner, but without modifying their positions they arent in cut off distance

    do counter3 = 1, processTracker
      transferThese(counter3, 1:3) = transferThese(counter3, 1:3) - lowerAndUpper(1:3)
    end do 

    if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3),0.0D+0]
    if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3),0.0D+0]
    if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6),0.0D+0]

    do counter3 = 1, oppositeCount
      transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) - currParticle(1:3)
    end do 

    !Sends/Receives the particles to the processors
    CALL MPI_Sendrecv(transferThese(1:processTracker, 1:4), processTracker*4, MPI_DOUBLE, adjacents(counter1), 10, &
                    receivedThese(1:numRecieved, 1:4), numRecieved*4, MPI_DOUBLE, adjacents(counter1+1), 10, &
                MPI_COMM_WORLD, status1, ierr)

    call MPI_Sendrecv(transferTheseAfter(1:oppositeCount, 1:4), oppositeCount*4, MPI_DOUBLE,  adjacents(counter1+1), 10, &
            receivedTheseAfter(1:numRecievedAfter, 1:4), numRecievedAfter*4, MPI_DOUBLE,  adjacents(counter1), 10, &
        & MPI_COMM_WORLD, status1, ierr)

    !Adds the particle offsets received to the relevant boundary vertex to get the position of the particle in local coordinate space
    if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3),0.0D+0]
    if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3),0.0D+0]
    if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6),0.0D+0]

    do counter3 = 1, numRecieved
      receivedThese(counter3, 1:3) = receivedThese(counter3, 1:3) + currParticle(1:3)
    end do

    do counter3 = 1, numRecievedAfter
      receivedTheseAfter(counter3, 1:3) = receivedTheseAfter(counter3, 1:3) + lowerAndUpper(1:3)
    end do

    !Undoing the offset calculation for its particles it already had so the positions aren't wrong for future iterations
    do counter3 = 1, processTracker
      transferThese(counter3, 1:3) = transferThese(counter3, 1:3) + lowerAndUpper(1:3)
    end do 

    if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3),0.0D+0]
    if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3),0.0D+0]
    if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6),0.0D+0]

    do counter3 = 1, oppositeCount
      transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) + currParticle(1:3)
    end do 

    !if (rank == 0 .and. logTime == "Y") print "(a18,i1,a17,f20.16)", " [TIME] Iteration ", counter2, &
    !& " Communications: ", elapsed_time()
    times(5 + 4*(counter2-1) + 2) = elapsed_time()
    !Comparing the original particles it had to any newly received ones, but only in one direction to prevent double counting
    do i = 1, arraySize
      do j = 1, numRecieved
         particlePairCount = particlePairCount + inRange(sentProcessData(i, 1:3),&
         & receivedThese(j, 1:3), cutoff, rank)        
      end do
    end do

    !if (rank == 0 .and. logTime == "Y") print "(a18,i1,a16,f20.16)", " [TIME] Iteration ", counter2, &
    !& " Pair Counting: ", elapsed_time()
    times(5 + 4*(counter2-1) + 3) = elapsed_time()
    !Transfering all particle data into allDataReceived, then deallocating any RAM-Style arrays to be used in the future iteration.
    deallocate(allDataReceived)
    allocate(allDataReceived(processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter, 4))

    if (processTracker /= 0) allDataReceived(1:processTracker, 1:4) = transferThese(1:processTracker, 1:4)

    if (oppositeCount /= 0) allDataReceived(processTracker+1 : processTracker+oppositeCount, 1:4) &
    & = transferTheseAfter(1:oppositeCount, 1:4)

    if (leftoverTracker /= 0) allDataReceived(processTracker+oppositeCount+1 : processTracker+oppositeCount+leftoverTracker, 1:4) &
    & = addAtEnd(1:leftoverTracker, 1:4)

    if (numRecieved /= 0) allDataReceived(&
      &processTracker+oppositeCount+leftoverTracker+1: &
      & (processTracker+oppositeCount+leftoverTracker+numRecieved), &
      1:4) = receivedThese(1:numRecieved, 1:4)

    if (numRecievedAfter /= 0) allDataReceived(&
      &processTracker+oppositeCount+leftoverTracker+numRecieved+1: &
      & (processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter), &
      1:4) = receivedTheseAfter(1:numRecievedAfter, 1:4)


    deallocate(transferThese) 
    deallocate(transferTheseAfter) 
    deallocate(receivedThese) 
    deallocate(receivedTheseAfter) 
    deallocate(addAtEnd)

    !if (rank == 0 .and. logTime == "Y") print "(a18,i1,a18,f20.16,/)", " [TIME] Iteration ", counter2, &
    !& " Resizing Arrays: ", elapsed_time()
    times(5 + 4*(counter2-1) + 4) = elapsed_time()
    43 format(" [OUTP] Iteration ", i1, " complete.")
    !44 format(" [Time] Iteration ", i1, " completion:", f15.13)
    if (dologging .and. rank == 0) print 43, counter2
  end do 

  deallocate(sentProcessData)
  deallocate(allDataReceived)

  !Send particle data back to the main process for final counting

  call MPI_REDUCE(particlePairCount, rootCount, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! if (rank /= 0) then
  !   ! call MPI_SEND(particlePairCount, 1, MPI_INTEGER8, 0, 3, MPI_COMM_WORLD, ierr)
  ! else
  !   do i = 1, nprocs-1
  !     counter1 = 0
  !     call MPI_RECV(counter1, 1, MPI_INTEGER8, i, 3, MPI_COMM_WORLD, status1, ierr)
  !     particlePairCount = particlePairCount + counter1 
  !   end do
  if (rank == 0) print "(a47, f6.2, a3, i15, a7, //)", &
    & " [OUTP] Total number of unique pairs with cutoff ", cutoff, " : ", rootCount, " [OUTP]"
  !end if

  !if (rank == 0 .and. logTime == "Y") print *, "[Time] Main Program time:", timestamp()
  if (rank == 0) times(18) = timestamp()

  call MPI_REDUCE(times, rootTimes, 18, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
      print *, "[TIME] Average Times: "
      print *, "[TIME] Particle Generation: ", rootTimes(1)
      print *, "[TIME] Divide Axis: ", rootTimes(2)
      print *, "[TIME] Particle Distribution: ", rootTimes(3)
      print *, "[TIME] Getting Neighbour Processes: ", rootTimes(4)/nprocs
      print *, "[TIME] Pair Counting Own Particles: ", rootTimes(5)/nprocs
      print *, "[TIME] Iteration 1 Particle Region Calculation: ", rootTimes(6)/nprocs
      print *, "[TIME] Iteration 1 Communications: ", rootTimes(7)/nprocs
      print *, "[TIME] Iteration 1 Pair Counting: ", rootTimes(8)/nprocs
      print "(a37, f20.16, /)", " [TIME] Iteration 1 Resizing Arrays: ", rootTimes(9)/nprocs
      print *, "[TIME] Iteration 2 Particle Region Calculation: ", rootTimes(10)/nprocs
      print *, "[TIME] Iteration 2 Communications: ", rootTimes(11)/nprocs
      print *, "[TIME] Iteration 2 Pair Counting: ", rootTimes(12)/nprocs
      print "(a37, f20.16, /)", " [TIME] Iteration 2 Resizing Arrays: ", rootTimes(13)/nprocs
      print *, "[TIME] Iteration 3 Particle Region Calculation: ", rootTimes(14)/nprocs
      print *, "[TIME] Iteration 3 Communications: ", rootTimes(15)/nprocs
      print *, "[TIME] Iteration 3 Pair Counting: ", rootTimes(16)/nprocs
      print "(a37, f20.16, /)", " [TIME] Iteration 3 Resizing Arrays: ", rootTimes(17)/nprocs
      print *, "[TIME] Main Program Time: ", rootTimes(18)
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

    if (vectorMod(temp) < (cutoff*cutoff)) count = 1
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
  real(kind=8) function elapsed_time(get_total) result(elapsed)
    implicit none

    ! whether or not to return the running total
    logical, intent(in), optional :: get_total
    
    ! implicit static keeps values between function calls
    real(kind=8) :: start_time = 0.0, end_time = 0.0, running_total = 0.0
    logical :: initialised = .false.
  
    if (present(get_total) .and. initialised) then ! if asking for the total then return it
      elapsed = running_total
    else if (initialised) then ! if the function has already been called once before
      start_time = end_time
      end_time = MPI_WTIME()
      elapsed = end_time - start_time
      running_total = running_total + elapsed
    else ! if this is the first time the function is ran
      start_time = MPI_WTIME()
      initialised = .true.
      elapsed = 0.0
    end if
  end function elapsed_time

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
