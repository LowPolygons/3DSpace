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

    CALL random_seed(put=seed)
    CALL random_number(unit_particles)
    
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
  integer(kind=int64) :: rootCount, testCount
  integer :: processTracker, totalReceivedInt
  logical :: success, dologging
  character :: generateNewFiles, logging, writeToFiles, logTime
  integer :: particleProcRatio, arraySize, oppositeCount, leftoverTracker

  !Output-Relevant, variable sized arrays
  double precision, dimension(:,:), allocatable :: particlePositions                    !all particle positions
  double precision, dimension(:,:), allocatable :: processorRange, sentProcessData, processDataHolder
  double precision, dimension(:,:), allocatable :: transferThese, receivedThese, totalReceivedTemp, totalReceivedReal
  double precision, dimension(:,:), allocatable :: transferTheseAfter, receivedTheseAfter, allDataReceived 
  double precision, dimension(:,:), allocatable :: addAtEnd
  real(kind=8), dimension(4) :: times, rootTimes
  integer, dimension(8) :: seed
  !mpi stuff
  integer :: ierr, rank, nprocs, numParticles, numRecieved, numRecievedAfter
  integer, dimension(MPI_STATUS_SIZE) :: status1

  !link cell stuff
  double precision, dimension(3) :: cellSize
  integer, dimension(3) :: numCells, currentParticleCell, currCell, compareCell
  integer, dimension(:,:,:,:), allocatable :: cellHeads
  integer :: currentIndex, theirIndex
  logical :: allocated

  logical :: testing = .true. !If true, it skips the USER input section 

  !Init for MPI
  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)  

  if (rank == 0) then
    i = elapsed_time()
    j = timestamp()
    print *, "[Time] Initialised Timers"
  end if
  
  ! Confirming if new data should be created or not. Checks if relevant files exist, if not program ends with message
  CALL confirm_config_modify(rank, testing)

  i = timestamp()
  i = elapsed_time()
  times = 0.00000000000000d0

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL read_config(numParticles, cutoff, lowerBound, upperBound, logging, i, writeToFiles, logTime) 

  
  if (rank == 0) then
    !displays the configurations for the current run
    CALL display_config(dologging, i)

    !creates the particles based on the config and stores in particlePositions
    CALL create_particles(particlePositions, numParticles, lowerBound, upperBound, writeToFiles)

    !get the number of domains
    numDomains = determineNoSplits(nprocs)

    !change the cut off incase the boundaries require as such
    cutoff = min( min(cutoff, (upperBound(1)-lowerBound(1))/ (numDomains(1)*2)), &
      min( (upperBound(2)-lowerBound(2))/(numDomains(2)*2), (upperBound(3)-lowerBound(3))/(numDomains(3)*2)) )

    allocate(processorRange(nprocs,6))
    processorRange = getBounds(nprocs, lowerBound, upperBound, numDomains)
    
    if (dologging .and. rank == 0) print *, "[OUTP] Distributing particles among processors"

    !send the lower and upper domain to other processors
    do counter2 = 1, nprocs
      curr2(1:6) = processorRange(counter2, 1:6)
      if (counter2 == 1) then
        !Its upper and lower boundary
        lowerAndUpper(1:6) = curr2(1:6)
      else
        !boundary
        CALL MPI_SEND(curr2(1:6), 6, MPI_DOUBLE, counter2-1, (counter2-1)*2, MPI_COMM_WORLD, ierr)
      end if
    end do

  end if 

  if (rank /= 0) then
    CALL MPI_RECV(lowerAndUpper, 6, MPI_DOUBLE, 0, rank*2, MPI_COMM_WORLD, status1, ierr)
  end if

  !Getting important values to all processors
  CALL MPI_BCAST(cutoff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(numParticles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(numDomains, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  !get all particles into each process
  if (rank /= 0) allocate(particlePositions(numParticles, 4))
  CALL MPI_BCAST(particlePositions, numParticles*4, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

  !get the particles you care about
  CALL sort_own_particles(sentProcessData, particlePositions, lowerAndUpper, arraySize)

  times(1) = elapsed_time()
  if (dologging .and. rank == 0) print *, "[OUTP] Distribution Complete"

  adjacents = getAdjacents(rank, numDomains)

  !ensure all processors are at this point
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  allocate(allDataReceived(arraySize, 4))
  allDataReceived = sentProcessData

  if (dologging .and. rank == 0) print *, "[OUTP] Starting Particle Communications." 

  CALL PARTICLE_COMMUNICATIONS(allDataReceived)

  times(2) = elapsed_time()

  !gets the link cell data
  CALL ESTABLISH_CELL_DATA(numCells, cellSize, lowerAndUpper, cutOff, cellHeads, rank, dologging)

  !updates the indexes of allDataReceived and the cell header array
  CALL ESTABLISH_INDEXES(numCells, cellSize, lowerAndUpper, cutoff, cellHeads, allDataReceived, rank, dologging)

  !Now that all of the indexes and cellHeads have been sorted you need to actually count the pairs
  CALL CALCULATE_CELL_PAIRS(numCells, cellSize, cutoff, cellHeads, allDataReceived, rank, dologging, testCount)

  deallocate(sentProcessData)
  deallocate(allDataReceived)

  times(3) = elapsed_time()
  if (rank == 0) print *, "Link Cell pair counting finished"

  CALL MPI_REDUCE(testCount, rootCount, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_REDUCE(times, rootTimes, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
    if (logTime == "Y")  then
      print *, "[TIME] Total Program Time:", timestamp()
      print *, "[TIME] Average Particle Distribution:", rootTimes(1)/nprocs
      print *, "[TIME] Average Communications:", rootTimes(2)/nprocs
      print *, "[TIME] Link Cell Counting:", rootTimes(3)/nprocs
    end if
  end if 

  if (rank == 0) print "(a47, f6.2, a3, i15, a7, //)", &
    & " [OUTP] Total number of unique pairs with cutoff ", cutoff, " : ", rootCount, " [OUTP]"
  CALL MPI_FINALIZE(ierr)

contains
  subroutine calculate_cell_pairs(numCells, cellSize, cutoff, cellHeads, allDataReceived, rank, dologging, testCount)
    double precision, dimension(3), intent(in) :: cellSize
    integer, dimension(3), intent(in) :: numCells
    integer, dimension(3) :: currentParticleCell
    integer, dimension(:,:,:,:), allocatable, intent(in) :: cellHeads
    double precision, dimension(4) :: currParticle
    double precision, intent(in) :: cutoff
    double precision, dimension(:,:), intent(in) :: allDataReceived
    integer :: counter1, counter2, counter3, currentIndex, theirIndex, i
    logical :: allocated
    integer, intent(in) :: rank
    logical, intent(in) :: dologging
    !each link cell should add thie vector described by these and compare the pairs in that list
    integer, dimension(13) :: linkCellChecksX, linkCellChecksY, linkCellChecksZ
    integer(kind=int64), intent(out) :: testCount

    linkCellChecksX = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    linkCellChecksY = [-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]
    linkCellChecksZ = [-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1, 1]

    !Count the pairs with all relevant received particles
    if (dologging .and. rank == 0) print *, "[OUTP] Linked Cell Pair counting starting."  
    testCount = 0

    !only checking the ones it actually owns
    do counter1 = 2, numCells(1)+1
      do counter2 = 2, numCells(2)+1
        do counter3 = 2, numCells(3)+1
          ![counter1, counter2, counter3] is the current Cell
          currCell = [counter1, counter2, counter3]

          theirIndex = cellHeads(currCell(1),currCell(2),currCell(3),1)

          do i = 1, 13
            compareCell = get_correct_cell(currCell,[linkCellChecksX(i),linkCellChecksY(i),linkCellChecksZ(i)], numCells)
            if (compareCell(1) <= numCells(1)+2 .and. compareCell(2) <= numCells(2)+2 .and. compareCell(3) <= numCells(3)+2) then
            !now the compareCell must be compared against the curr Cell
            currentIndex = cellHeads(currCell(1), currCell(2), currCell(3), 1)
            theirIndex = cellHeads(compareCell(1), compareCell(2), compareCell(3), 1)

              do while(currentIndex /= -1)
                do while(theirIndex /= -1)
                  if (currentIndex /= theirIndex) then
                    testCount = testCount + inRange(allDataReceived(currentIndex, 1:4),&
                            & allDataReceived(theirIndex, 1:4), cutoff, rank)      
                  end if  
                  theirIndex = int(allDataReceived(theirIndex, 4)) 
                end do
                currentIndex = int(allDataReceived(currentIndex, 4))   
                theirIndex = cellHeads(compareCell(1), compareCell(2), compareCell(3), 1)        
              end do
            end if
          end do

          !itself
          currentIndex = cellHeads(currCell(1), currCell(2), currCell(3), 1)
          if (currentIndex /= - 1) theirIndex = int(allDataReceived(currentIndex, 4))
          
          !itself
          do while(currentIndex /= -1)
            do while(theirIndex /= -1)
              testCount = testCount + inRange(allDataReceived(currentIndex, 1:4),&
                      & allDataReceived(theirIndex, 1:4), cutoff, rank)       
              theirIndex = int(allDataReceived(theirIndex, 4)) 
            end do
            currentIndex = int(allDataReceived(currentIndex, 4))   
            if (currentIndex /= -1) theirIndex = int(allDataReceived(currentIndex, 4))
          end do
          !if (rank == 0) print *, ""
        end do
      end do 
    end do
  end subroutine calculate_cell_pairs

  subroutine establish_indexes(numCells, cellSize, lowerAndUpper, cutoff, cellHeads, allDataReceived, rank, dologging)
    double precision, dimension(3), intent(in) :: cellSize
    integer, dimension(3), intent(in) :: numCells
    integer, dimension(3) :: currentParticleCell
    integer, dimension(:,:,:,:), allocatable, intent(inout) :: cellHeads
    double precision, dimension(6) :: lowerAndUpper
    double precision, dimension(4) :: currParticle
    double precision, intent(in) :: cutoff
    double precision, dimension(:,:), intent(inout) :: allDataReceived
    integer :: counter1, currentIndex, theirIndex, i
    logical :: allocated
    integer, intent(in) :: rank
    logical, intent(in) :: dologging

    !THis array will update the indexes in the all particle positions loop to maintain consistent indexes
    do counter1 = 1, size(allDataReceived)/4
      !first, do the offset, as well as completely localising it to be around 0
      allDataReceived(counter1, 1) = allDataReceived(counter1, 1) - lowerAndUpper(1) + cellSize(1) !offsetting all the particles by the cellSize to ensure they all fall within the correct boundary
      allDataReceived(counter1, 2) = allDataReceived(counter1, 2) - lowerAndUpper(2) + cellSize(2) !offsetting all the particles by the cellSize to ensure they all fall within the correct boundary
      allDataReceived(counter1, 3) = allDataReceived(counter1, 3) - lowerAndUpper(3) + cellSize(3) !offsetting all the particles by the cellSize to ensure they all fall within the correct boundary

      currParticle(1:3) = allDataReceived(counter1, 1:3)

      currentIndex = counter1 !they need new indexes

      !their indexes will be all over the place so this localises them initally
      allDataReceived(counter1, 4) = currentIndex !just in case

      !current cell it is in, +1, because the indexes start at 1
      currentParticleCell(1) = INT((currParticle(1))/cellSize(1)) + 1
      currentParticleCell(2) = INT((currParticle(2))/cellSize(2)) + 1
      currentParticleCell(3) = INT((currParticle(3))/cellSize(3)) + 1

      !All code in the loop from here checks if there are any items in the list, and if there are it goes along the path
      !otherwise, it sets the new head particle index

      if (cellHeads(currentParticleCell(1), currentParticleCell(2), currentParticleCell(3), 1) == -1) then
        cellHeads(currentParticleCell(1), currentParticleCell(2), currentParticleCell(3), 1) = currentIndex
        allDataReceived(counter1, 4) = DBLE(-1)
        allocated = .true.
      else
        theirIndex = cellHeads(currentParticleCell(1), currentParticleCell(2), currentParticleCell(3), 1)
        allocated = .false.
      end if 

      do while (.not. allocated)
        !therefore it needs to be allocated, theirIndex already has the index of the head particle
        i = int(allDataReceived(theirIndex, 4))

        if (i == -1) then
          !reached the end of that part of the list
          allocated = .true.
          allDataReceived(theirIndex, 4) = currentIndex
          allDataReceived(counter1, 4) = -1
        else
          theirIndex = i
        end if 
      end do

      allDataReceived(counter1, 1) = allDataReceived(counter1, 1) + lowerAndUpper(1) - cellSize(1)
      allDataReceived(counter1, 2) = allDataReceived(counter1, 2) + lowerAndUpper(2) - cellSize(2)
      allDataReceived(counter1, 3) = allDataReceived(counter1, 3) + lowerAndUpper(3) - cellSize(3)
    end do
  end subroutine establish_indexes

  subroutine establish_cell_data(numCells, cellSize, lowerAndUpper, cutoff, cellHeads, rank, dologging)
    double precision, dimension(3), intent(out) :: cellSize
    integer, dimension(3), intent(out) :: numCells
    integer, dimension(3) :: currentParticleCell
    integer, dimension(:,:,:,:), allocatable, intent(out) :: cellHeads
    double precision, dimension(6) :: lowerAndUpper
    double precision, intent(in) :: cutoff

    integer, intent(in) :: rank
    logical, intent(in) :: dologging

    numCells(1) = INT((lowerAndUpper(4)-lowerAndUpper(1))/(cutoff)) 
    numCells(2) = INT((lowerAndUpper(5)-lowerAndUpper(2))/(cutoff))
    numCells(3) = INT((lowerAndUpper(6)-lowerAndUpper(3))/(cutoff))

    !note: this is the number of cells without including neighbours
    cellSize(1) = (lowerAndUpper(4)-lowerAndUpper(1)) / DBLE( numCells(1) ) !casts to and int and back tto a double to automatiCALLy truncate
    cellSize(2) = (lowerAndUpper(5)-lowerAndUpper(2)) / DBLE( numCells(2) ) !casts to and int and back tto a double to automatiCALLy truncate
    cellSize(3) = (lowerAndUpper(6)-lowerAndUpper(3)) / DBLE( numCells(3) ) !casts to and int and back tto a double to automatiCALLy truncate


    if (rank == 0 .and. dologging) then
      print *, "[OUTP] Difference: ", lowerAndUpper(4:6)-lowerAndUpper(1:3)
      print *, "[OUTP] Num Cells: ", numCells
      print *, "[OUTP] Cell Size: ", cellSize 
      print *, "[OUTP] Does it fit: ", numCells*cellSize
    end if 

    allocate( cellHeads( numCells(1)+2, numCells(2)+2, numCells(3)+2, 1))

    cellHeads = -1 !-1 is empty
  end subroutine establish_cell_data
  subroutine particle_communications(allDataReceived)
    implicit none
    double precision, dimension(:,:), allocatable, intent(inout) :: allDataReceived
    
    integer :: processTracker, oppositeCount, leftoverTracker
    double precision, dimension(:,:), allocatable :: transferThese, receivedThese
    double precision, dimension(:,:), allocatable :: transferTheseAfter, receivedTheseAfter
    double precision, dimension(:,:), allocatable :: addAtEnd
    double precision, dimension(4) :: currParticle
    integer :: counter1, counter2
    !The current iteration
    
    counter2 = 0
    
    do counter1 = 1, 5, 2
      
      counter2 = counter2 + 1 
      allocate(transferThese(size(allDataReceived)/4, 4))
      allocate(transferTheseAfter(size(allDataReceived)/4, 4)) 
      allocate(addAtEnd(size(allDataReceived)/4, 4))

      processTracker = 0 
      oppositeCount = 0 
      leftoverTracker = 0 
      numRecieved = 0
      numRecievedAfter = 0

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
      
      CALL MPI_Sendrecv(processTracker, 1, MPI_INTEGER, adjacents(counter1), 2, &
                          numRecieved, 1, MPI_INTEGER, adjacents(counter1+1), 2, &
                  MPI_COMM_WORLD, status1, ierr)


      CALL MPI_Sendrecv(oppositeCount, 1, MPI_INTEGER, adjacents(counter1+1), 3, &
                      numRecievedAfter, 1, MPI_INTEGER, adjacents(counter1), 3, &
                    MPI_COMM_WORLD, status1, ierr)

      allocate(receivedThese(numRecieved, 4))
      allocate(receivedTheseAfter(numRecievedAfter, 4))


      do counter3 = 1, processTracker
        transferThese(counter3, 1:3) = transferThese(counter3, 1:3) - lowerAndUpper(1:3)
      end do 

      if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3),0.0D+0]
      if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3),0.0D+0]
      if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6),0.0D+0]

      do counter3 = 1, oppositeCount
        transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) - currParticle(1:3)
      end do 

      CALL MPI_Sendrecv(transferThese(1:processTracker, 1:4), processTracker*4, MPI_DOUBLE, adjacents(counter1), 10, &
                      receivedThese(1:numRecieved, 1:4), numRecieved*4, MPI_DOUBLE, adjacents(counter1+1), 10, &
                  MPI_COMM_WORLD, status1, ierr)

      CALL MPI_Sendrecv(transferTheseAfter(1:oppositeCount, 1:4), oppositeCount*4, MPI_DOUBLE,  adjacents(counter1+1), 10, &
              receivedTheseAfter(1:numRecievedAfter, 1:4), numRecievedAfter*4, MPI_DOUBLE,  adjacents(counter1), 10, &
          & MPI_COMM_WORLD, status1, ierr)

      if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3),0.0D+0]
      if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3),0.0D+0]
      if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6),0.0D+0]

      do counter3 = 1, numRecieved
        receivedThese(counter3, 1:3) = receivedThese(counter3, 1:3) + currParticle(1:3)
      end do

      do counter3 = 1, numRecievedAfter
        receivedTheseAfter(counter3, 1:3) = receivedTheseAfter(counter3, 1:3) + lowerAndUpper(1:3)
      end do

      do counter3 = 1, processTracker
        transferThese(counter3, 1:3) = transferThese(counter3, 1:3) + lowerAndUpper(1:3)
      end do 

      if (counter2 == 1) currParticle = [lowerAndUpper(4),lowerAndUpper(2),lowerAndUpper(3),0.0D+0]
      if (counter2 == 2) currParticle = [lowerAndUpper(1),lowerAndUpper(5),lowerAndUpper(3),0.0D+0]
      if (counter2 == 3) currParticle = [lowerAndUpper(1),lowerAndUpper(2),lowerAndUpper(6),0.0D+0]

      do counter3 = 1, oppositeCount
        transferTheseAfter(counter3, 1:3) = transferTheseAfter(counter3, 1:3) + currParticle(1:3)
      end do 

      deallocate(allDataReceived)
      allocate(allDataReceived(processTracker+oppositeCount+leftoverTracker+numRecieved+numRecievedAfter, 4))

      if (processTracker /= 0) allDataReceived(1:processTracker, 1:4) = transferThese(1:processTracker, 1:4)

      if (oppositeCount /= 0) allDataReceived(processTracker+1 : processTracker+oppositeCount, 1:4) &
      & = transferTheseAfter(1:oppositeCount, 1:4)

      if (leftoverTracker /= 0) &
      & allDataReceived(processTracker+oppositeCount+1 : processTracker+oppositeCount+leftoverTracker, 1:4) &
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
    end do
  end subroutine particle_communications

  subroutine sort_own_particles(sentProcessData, particlePositions, lowerAndUpper, arraySize)
    double precision, dimension(:,:), allocatable, intent(out) :: sentProcessData
    double precision, dimension(:,:), allocatable :: processDataHolder
    double precision, dimension(:,:), intent(in) :: particlePositions
    double precision, dimension(6), intent(in) :: lowerAndUpper
    integer :: processTracker
    integer, intent(out) :: arraySize
    double precision, dimension(4) :: currParticle

    allocate(processDataHolder(size(particlePositions)/4, 4))

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
  end subroutine sort_own_particles

  function get_correct_cell(currentCell, sumVector, numCellsTotal) result(newCell)
    integer, dimension(3), intent(in) :: currentCell, sumVector, numCellsTotal
    integer, dimension(3) :: newCell
    
    newCell = currentCell + sumVector
    !x
    if (newCell(1) > (numCellsTotal(1)+2)) newCell(1) = (numCellsTotal(1)+2)-newCell(1)
    if (newCell(1) < 1) newCell(1) = numCellsTotal(1)+2
    !y
    if (newCell(2) > (numCellsTotal(2)+2)) newCell(2) = (numCellsTotal(2)+2)-newCell(2)
    if (newCell(2) < 1) newCell(2) = numCellsTotal(2)+2
    !z
    if (newCell(3) > (numCellsTotal(3)+2)) newCell(3) = (numCellsTotal(3)+2)-newCell(3)
    if (newCell(3) < 1) newCell(3) = numCellsTotal(3)+2
  end function get_correct_cell
  !Simply gets the modulus of the distance and compares against the cutoff, nothing fancy like PBC
  subroutine create_particles(particlePositions, numParticles, lowerBound, upperBound, writeToFiles)

    double precision, dimension(:,:), intent(out), allocatable :: particlePositions
    integer, intent(in) :: numParticles
    double precision, dimension(3), intent(in) :: lowerBound, upperBound
    character, intent(in) :: writeToFiles
    logical :: success

    ! All particle positions, only known for this process at this time
    allocate(particlePositions(numParticles, 4))
    
    inquire(file="coordinates.txt", exist=success)
    if (success .and. writeToFiles == "Y") then
      open(11, file="coordinates.txt", status="old")

      do counter1 = 1, numParticles
        read(11,*) particlePositions(counter1, 1), particlePositions(counter1, 2), particlePositions(counter1, 3)
        particlePositions(counter1, 4) = counter1
      end do 

      close(11)
    else 
      CALL random_number(particlePositions)
      !Formatting the particles to be distributed evenly across the given upper and lower bound
      do counter1 = 1, numParticles
        particlePositions(counter1, 1:3) = [&
          & lowerBound(1) + particlePositions(counter1, 1)*(upperBound(1) - lowerBound(1)), &
          & lowerBound(2) + particlePositions(counter1, 2)*(upperBound(2) - lowerBound(2)), &
          & lowerBound(3) + particlePositions(counter1, 3)*(upperBound(3) - lowerBound(3))  &
        &]
      end do
    end if 
  end subroutine create_particles

  subroutine confirm_config_modify(rank, testing)
    integer :: rank
    logical :: testing, success

    if (rank == 0 .and. .not. testing) then
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
  end subroutine confirm_config_modify

  subroutine display_config(dologging, seedVal)
    logical, intent(out) :: dologging
    integer, dimension(8) :: seed
    integer :: seedVal
    
    seed = seedVal

    CALL random_seed(put=seed)

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
  end subroutine display_config

  subroutine read_config(numParticles, cutoff, lowerBound, upperBound, logging, seed, writeToFiles, logTime) 
    integer, intent(out) :: numParticles, seed
    double precision, intent(out) :: cutoff
    double precision, dimension(3), intent(out) :: lowerBound, upperBound
    character, intent(out) :: logging, writeToFiles, logTime

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
  end subroutine read_config
  
  integer function inRange(p1, p2, cutoff, rank) result(count)
    double precision, dimension(1:4), intent(in) :: p1, p2
    double precision, dimension(1:4) :: temp
    double precision, intent(in):: cutoff
    integer, intent(in) :: rank

    temp = p2 - p1
    count = 0

    if (vectorMod(temp) < (cutoff*cutoff)) count = 1

    !if (count == 1) print *, rank, p1, p2
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
    
    ! implicit static keeps values between function CALLs
    real(kind=8) :: start_time = 0.0, end_time = 0.0, running_total = 0.0
    logical :: initialised = .false.
  
    if (present(get_total) .and. initialised) then ! if asking for the total then return it
      elapsed = running_total
    else if (initialised) then ! if the function has already been CALLed once before
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
