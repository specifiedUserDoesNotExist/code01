!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! PRESS, William H.; TEUKOLSKY, Saul A. Numerical calculation of derivatives. 
! Computers in Physics, v. 5, n. 1, p. 68-69, 1991.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine derivada(x, y, p, s, ans, erro)
	! s = 1 
	!	derivada dx2
	! s = 2
	!	derivada dy2
	! s = 3
	!	derivada dxdy
	! p = chi2fun(x,y)

	use global
	
	implicit none
	integer, parameter :: N = 15
	integer :: s
	real (kind = 8) :: x, y, p
	real (kind = 8):: erro, errAux, ans, CON2, delta, chi2fun
	real (kind = 8), parameter :: CON = 1.4D0
	real (kind = 8), parameter :: SAFE = 2.D0
	
	real (kind = 8) :: a1, a2
	
	integer :: i, j
	real (kind = 8) :: vec(N, N)
	
	erro = 1000.d0
	delta = 0.001d0
	
	!wm = x
	!wl = y
	
	!print*, chi2fun(wm, wl)
	
	if (s == 1) then
		wm = x + delta
		wl = y
		vec(1,1) = chi2fun(wm, wl)

		wm = x - delta
		vec(1,1) = (vec(1,1) + chi2fun(wm, wl) - 2*p)/(delta*delta)
		
		!vec(1,1) = (chi2fun(wm, wl) + chi2fun(wm-delta, wl) - 2*p)/(delta*delta)
	else if (s == 2) then
		wm = x
		wl = y + delta
		vec(1,1) = chi2fun(wm, wl)
		
		wl = y - delta
		vec(1,1) = (vec(1,1) + chi2fun(wm, wl) - 2*p)/(delta*delta)
		
		!vec(1,1) = (chi2fun(wm, wl+delta) + chi2fun(wm, wl-delta) - 2*p)/(delta*delta)
	else
		wm = x + delta
		wl = y + delta
		vec(1,1) = chi2fun(wm, wl)
		
		wl = y - delta
		vec(1,1) = vec(1,1) - chi2fun(wm, wl)
		
		wm = x - delta
		vec(1,1) = vec(1,1) + chi2fun(wm, wl)
		
		wl = y + delta
		vec(1,1) = (vec(1,1) - chi2fun(wm, wl))/(4*delta*delta)
		
		
		!vec(1,1) = ((chi2fun(wm+delta, wl+delta)+chi2fun(wm-delta, wl-delta)) - &
		!(chi2fun(wm+delta, wl-delta)+chi2fun(wm-delta, wl+delta))) &
		!/(4*delta*delta)
	end if

	do i = 2 , N
		delta = delta/CON
		!vec(1, i) = (chi2fun(wm+delta) + chi2fun(wm-delta) - 2*chi2fun(wm))/(delta*delta)
		if (s == 1) then
			wm = x + delta
			wl = y
			vec(1,i) = chi2fun(wm, wl)

			wm = x - delta
			vec(1,i) = (vec(1,i) + chi2fun(wm, wl) - 2*p)/(delta*delta)
			
			!vec(1,i) = (chi2fun(wm, wl) + chi2fun(wm-delta, wl) - 2*p)/(delta*delta)
		else if (s == 2) then
			wm = x
			wl = y + delta
			vec(1,i) = chi2fun(wm, wl)
			
			wl = y - delta
			vec(1,i) = (vec(1,i) + chi2fun(wm, wl) - 2*p)/(delta*delta)
			
			!vec(1,1) = (chi2fun(wm, wl+delta) + chi2fun(wm, wl-delta) - 2*p)/(delta*delta)
		else
			wm = x + delta
			wl = y + delta
			vec(1,i) = chi2fun(wm, wl)
			
			wl = y - delta
			vec(1,i) = vec(1,i) - chi2fun(wm, wl)
			
			wm = x - delta
			vec(1,i) = vec(1,i) + chi2fun(wm, wl)
			
			wl = y + delta
			vec(1,i) = (vec(1,i) - chi2fun(wm, wl))/(4*delta*delta)
			
			
			!vec(1,1) = ((chi2fun(wm+delta, wl+delta)+chi2fun(wm-delta, wl-delta)) - &
			!(chi2fun(wm+delta, wl-delta)+chi2fun(wm-delta, wl+delta))) &
			!/(4*delta*delta)
		end if

		CON2 = CON*CON
		do j = 2, i
			vec(j, i) = (vec(j-1, i)*CON2 - vec(j-1, i-1))/(CON2 - 1)
			CON2 = CON2 * CON2
			errAux = max(abs(vec(j, i) - vec(j-1, i)), abs(vec(j, i) - vec(j-1, i-1)))
			
			if (errAux <= erro) then
				erro = errAux
				ans = vec(j, i)
			end if
		end do
		
		if (abs(vec(i, i) - vec(i-1, i-1)) >= SAFE*erro) then
			return
		end if
	end do
	
	return
	
end subroutine
