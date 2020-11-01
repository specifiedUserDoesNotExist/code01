subroutine readNumberLines(fileIn, nlines)
	
	implicit none
	integer :: io
	integer, intent(out) :: nlines
	
	character (len = 30), intent(in) :: fileIn
	
	open(1, file=fileIn)
	nlines = 0
	do 
		read (1, *, iostat=io)
		if (io/=0) exit
		nlines = nlines + 1
	end do
	close(1)
	
end subroutine readNumberLines
