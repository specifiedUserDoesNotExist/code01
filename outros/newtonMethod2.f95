module global
	implicit none
	
	real (kind = 8) :: N, l
	real (kind = 8) :: k
	
end module

subroutine derivada(x, ans, erro)
	
	use global
	implicit none

	integer, parameter :: nn = 15
	real (kind = 8) :: x
	real (kind = 8) :: erro, errAux, ans, CON2, delta, func
	real (kind = 8), parameter :: CON = 1.400
	
	integer :: i,j
	real (kind = 8) :: vec(nn,nn)
	
	erro = 1000.d0
	delta = 0.001d0
	
	vec(1,1) = (func(x+delta)-func(x-delta))/(2.0*delta)
	
	
	do i=2,nn
		delta = delta/CON
		vec(1,i) = (func(x+delta)-func(x-delta))/(2.0*delta)
		
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

subroutine newtonsMethod(x0, x1)

	use global
	implicit none

	real (kind = 8) :: x0, x1, eps, erro
	real (kind = 8) :: f0, f1, func
	integer, parameter :: interations = 1000
	integer :: i
	
	eps = 1e-7
  
	do i=1, interations
		f0 = func(x0)
		call derivada(x0, f1, erro)
		
		if (abs(f1) < eps) then
			print*, "nao converge"
			return
		end if
		
		x1 = x0 - f0/f1
		
		if (abs(x1-x0) <= eps) then
			return
		end if
		
		x0 = x1
	
	end do

end subroutine

function func(x)

	use global
	implicit none
	
	real (kind = 8) :: x, func
	
	func = 2.0*l*N - x*sqrt(1.0 + (l*k*x)**2.) - x
	
end function

function ns(f)

	use global
	implicit none
	
	real (kind = 8) :: ns, f
	real (kind = 8) :: lkf, n1, d1, n2, d2
	
	lkf = l*k*f
	
	n1 = 1+2*k
	d1 = 1+lkf**2
	
	n2 = 2*k
	d2 = (d1**3)*(sqrt(d1)-lkf)

	ns = 1 - l**2*(n1/d1 - n2/d2)

end function

program main

	use global
	implicit none
	
	real (kind = 8) :: ns, f
	real (kind = 8) :: x = 0
	integer :: INT = 1000
	integer :: ki
	
	N = 60.0
	l = 1.5
	
	do ki = 1, INT
		k = ki/(1.d0*INT)
		
		call newtonsMethod(x, f)
		
		print*, k, f, ns(f)
	end do

end program
