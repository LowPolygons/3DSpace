program mpiTest
    use mpi

    integer :: rank, nprocs, ierr

    !status variable tells the status of send/receive calls, needed for receive subroutine
    integer, dimension(MPI_STATUS_SIZE) :: status1

    !arrangement for saving the host name
    character*(MPI_MAX_PROCESSOR_NAME) :: hostname
    integer :: namesize 

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

    !by default these were all the same
    print *, "Hello, I am ", hostname(1:namesize), " with rank ", rank, " of ", nprocs, "process"

    !SE
    call MPI_FINALIZE(ierr)

end program mpiTest

!  mpirun -n 10 ./test
!  Hello, I am DLNLAB0010 with rank            4  of           10 process
!  Hello, I am DLNLAB0010 with rank            5  of           10 process
!  Hello, I am DLNLAB0010 with rank            6  of           10 process
!  Hello, I am DLNLAB0010 with rank            7  of           10 process
!  Hello, I am DLNLAB0010 with rank            9  of           10 process
!  Hello, I am DLNLAB0010 with rank            0  of           10 process
!  Hello, I am DLNLAB0010 with rank            1  of           10 process
!  Hello, I am DLNLAB0010 with rank            2  of           10 process
!  Hello, I am DLNLAB0010 with rank            3  of           10 process
!  Hello, I am DLNLAB0010 with rank            8  of           10 process