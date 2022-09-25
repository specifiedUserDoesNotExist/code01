subroutine derivada(x, k, ans, erro)

	implicit none
	integer, parameter :: N = 15
	real (kind = 8) :: x, k
	real (kind = 8) :: erro, errAux, ans, CON2, delta, func
	real (kind = 8), parameter :: CON = 1.400
	
	integer :: i,j
	real (kind = 8) :: vec(N,N)
	
	erro = 1000.d0
	delta = 0.001d0
	
	vec(1,1) = (func(x+delta, k)-func(x-delta, k))/(2.0*delta)
	
	
	do i=2,N
		delta = delta/CON
		vec(1,i) = (func(x+delta, k)-func(x-delta, k))/(2.0*delta)
		
		CON2 = CON*CON
		
		do j = 2, i
			vec(j,i) = (vec(j-1,i)*CON2-vec(j-1,i-1))/(CON2-1.)
			CON2 = CON2*CON2
			errAux = max(abs(vec(j,i)-vec(j-1,i)),abs(vec(j,i)-vec(j-1,i-1)))
			
			if (errAux <= erro) then
				erro = errAux
				ans = vec(j,i)
			end if
		end do
		
		if (abs(vec(i,i)-vec(i-1,i-1)) >= 2*erro) then
			return
		end if
	end do
	
	return
end subroutine

subroutine newtonsMethod(x0, k, x1)

	implicit none
	real (kind = 8) :: x0, x1, eps, erro
	real (kind = 8) :: f0, f1, func, k
	integer, parameter :: interations = 1000
	integer :: i
	
	eps = 1e-7
	!x0=0.0
	
	do i=1, interations
		f0 = func(x0, k)
		call derivada(x0, k, f1, erro)
		
		if (abs(f1) < eps) then
			return
		end if
		
		x1 = x0 - f0/f1
		
		if (abs(x1-x0) <= eps) then
			return
		end if
		
		x0 = x1
	
	end do
		
	print*, 'Not converged'
	x1=0.0
	return

end subroutine

function func(x, k)

	implicit none
	
	real (kind = 8) :: x, func
	real (kind = 8) :: N, l, k
	
	N = 60.0
	l = 1.5
	
	!func = x**2 - 5*x + 6  
	func = 2.0*l*N - x*sqrt(1.0 + (l*k*x)**2.) - x
	
end function

program main
	implicit none
	
	real (kind = 8) :: x = 0
	real (kind = 8) ans, erro, kaux
	integer :: INT = 10000
	integer :: k
	
	do k =1, INT
		kaux = k/(1.d0*INT)
		
		call newtonsMethod(x,kaux, ans)
		
		print*, kaux, ans
	end do

end program
