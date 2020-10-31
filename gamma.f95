!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!          GAMMA DISTRIBUTION FUNCTION
!     (C) Copr. 1986-92 Numerical Recipes Software .U:!0)j3.
!-------------------------------------------------------------------------------
!------------------------ FUNCAO GAMA INCOMPLETA -------------------------------
!-------------------------------------------------------------------------------
function gammp(a, x)
	
	real (kind = 8) :: a, gammp, x
	real (kind = 8) :: gammcf, gamser, gln
	
	if (x < 0 .or. a <= 0) then
		write (*,*) "a = ", a
		write (*,*) "x = ", x
		print *, "argumentos ruins na gammp"
		call exit(1)
	end if
	
	if ( x < a+1.) then
		call gser(gamser, a, x, gln)
		gammp = gamser
	else
		call gcf(gammcf, a, x, gln)
		gammp = 1. - gammcf
	endif

end function gammp
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine gser(gamser, a, x, gln)

	integer, parameter :: ITMAX = 100
	integer :: n
	
	real (kind = 8) :: a, gamser, gln, x
	real (kind = 8), parameter :: EPS = 3.d-7
	real (kind = 8) :: ap, del, sum, gammln 

	gln = gammln(a)
         
	if (x <= 0.) then
		
		if (x < 0.) then
			print *, 'x < 0 no gser'
			call exit(1)
		end if
		
		gamser = 0.
		return 
	endif
         
	ap = a
	sum = 1./a
	del = sum
         
	do n = 1, ITMAX
		
		ap=ap+1.
		del = del*x/ap
		sum = sum+del
		
		if ( abs(del) < abs(sum)*EPS ) then
			gamser = sum*exp(-x+a*log(x)-gln)
			return
		end if 
	end do
	
	!if ( abs(del) > abs(sum)*EPS ) then
	print *, "a muito grande, ITMAX muito pequeno no gser"
	call exit(1)
	!end if
end subroutine gser
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine gcf(gammcf, a, x, gln)

	integer , parameter:: ITMAX = 100
	integer :: i
	
	real (kind = 8) :: a, gammcf, gln, x
	real (kind = 8) :: an, b, c, d, del, h, gammln
	real (kind = 8), parameter :: EPS = 3.d-7
	real (kind = 8), parameter :: FPMIN = 1.d-30
	
	gln = gammln(a)
	b = x+1. -a
	c = 1./FPMIN
	d = 1./b
	h = d

	do i=1, ITMAX
           
		an =-i*(i-a)
		b=b+2
		
		d =an*d+b
		if (abs(d) < FPMIN ) then 
			d = FPMIN
		end if
		
		c=b+an/c
		if(abs(c) < FPMIN) then 
			c = FPMIN
		end if
		
		d = 1./d
		del = d*c
		h = h*del
		
		if (abs(del-1.) < EPS) then
			gammcf= exp(-x+a*log(x)-gln)*h
			return
		end if
	end do
         
	!if ( abs(del-1.) > EPS ) then 
	print *, "a grande, ITMAX pequeno no gcf"
	call exit(1)
	!end if
end subroutine gcf
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
function gammln(xx)

	real (kind = 8) :: gammln, xx
	real (kind = 8) :: ser, tmp, x, y
	integer :: j
	
	real (kind = 8), save :: stp = 2.5066282746310005d0
	real (kind = 8), save :: cof(6) = (/76.18009172947146d0, &
			-86.50532032941677d0,	&
			24.01409824083091d0,	& 
			-1.231739572450155d0,	&
			.1208650973866179d-2,	&
			-.5395239384953d-5/)
	
	x = xx
	y = x
	tmp = x + 5.5d0
	tmp = (x + 0.5d0)*log(tmp)-tmp
	ser = 1.00000000090015d0
	
	do j=1, 6
		y=y+1.d0
		ser=ser+cof(j)/y
	end do
	
	gammln = tmp+log(stp*ser/x)
end function gammln
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!      (C) Copr. 1986-92 Numerical Recipes Software .U:!0)j3.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  FIM   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
