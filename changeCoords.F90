program coordsGen
	implicit none
	
	double precision, dimension(300000, 3) :: particlePositions !domain 000,111
	double precision, dimension(3) :: theOffset = [0.00000000000000d0, 0.00000000000000d0, 1.0000000000000d0]
	integer :: counter1 
	
    !open(11, file="coordinates.txt", status="old")
    open(11, file="new.txt", status="old")
    !open(11, file="originals.txt", status="old")


    do counter1 = 1, 300000
		read(11,*) particlePositions(counter1, 1), particlePositions(counter1, 2), particlePositions(counter1, 3)
    end do 

    close(11)
	
	
	open(12, file="new.txt", status="old")
	
	9 format(f20.17, 4x, f20.17, 4x, f20.17)

	do counter1 = 1, 300000
		write(12, 9) (particlePositions(counter1, 1) + theOffset(1)), (particlePositions(counter1, 2) + theOffset(2)), &
		&(particlePositions(counter1,3) + theOffset(3))
	end do 
	
	close(12)
	
	print *, "done"
end program coordsGen