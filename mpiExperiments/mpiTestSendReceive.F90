program mpiTest
    use mpi

    integer :: rank, nprocs, ierr

    !status variable tells the status of send/receive calls, needed for receive subroutine
    integer, dimension(MPI_STATUS_SIZE) :: status1

    !arrangement for saving the host name
    character*(MPI_MAX_PROCESSOR_NAME) :: hostname
    integer :: namesize, i, sumofRanks, temp
    
    sumofRanks = 0
    !SE
    call MPI_INIT(ierr)

    !set up communicator size
	!var order :: Communicator, NumProcessers, Error variable
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	!set up ranks
	!var order :: Communicator, Rank/ProcessoridVariable, Error variable
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)	

    !get the hostname for each process
    call MPI_GET_PROCESSOR_NAME(hostname, namesize, ierr)

    ! Starting master process

    !call random_number(addThese)
   
    print *, "Process: ", rank

    if ( rank /= 0 ) then
        !send data to other ranks
        !syntax :: call MPI_SEND(start_address/variablebeingsent, count/numDataPieces, datatype, destination pid/rank, tag, communicator, ierr)
        !Tag is a unique identifier, tag pairs up send and receive calls
        !integers are 4 bytes in MPI_INT

        call MPI_SEND(rank, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, ierr)
        print *, "Data sent from Process", rank, "to Master Process."
    else !non master process
        print *, nprocs

        do i = 1, nprocs-1
            !receive data from ranks
            !syntax: call MPI_RECV(startAddress, count, datatype, source pid/rank, tag, communicator, status, ierr)

            call MPI_RECV(temp, 1, MPI_INT, i, i, MPI_COMM_WORLD, status1, ierr)
            print *, "Iteration: ", i, ", Received :", temp
            sumofRanks = sumofRanks + temp
        end do

        print *, "Sum: ", sumofRanks
    end if

    !SE
    call MPI_FINALIZE(ierr)

end program mpiTest

!  Process:            0
!           10
!  Iteration:            1 , Received :           1
!  Iteration:            2 , Received :           2
!  Iteration:            3 , Received :           3
!  Iteration:            4 , Received :           4
!  Iteration:            5 , Received :           5
!  Iteration:            6 , Received :           6
!  Iteration:            7 , Received :           7
!  Iteration:            8 , Received :           8
!  Iteration:            9 , Received :           9
!  Sum:           45
!  Process:            1
!  Data sent from Process           1 to Master Process.
!  Process:            2
!  Data sent from Process           2 to Master Process.
!  Process:            4
!  Data sent from Process           4 to Master Process.
!  Process:            5
!  Data sent from Process           5 to Master Process.
!  Process:            6
!  Data sent from Process           6 to Master Process.
!  Process:            7
!  Data sent from Process           7 to Master Process.
!  Process:            8
!  Data sent from Process           8 to Master Process.
!  Process:            9
!  Data sent from Process           9 to Master Process.
!  Process:            3
!  Data sent from Process           3 to Master Process.