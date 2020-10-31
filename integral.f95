!-------------------------------------------------------------------------------
! ------------------------------FUNCAO EXTERNA----------------------------------
!-------------------------------------------------------------------------------
!--------------------------- ROTINA DE INTEGRACAO ------------------------------
subroutine qromb1 (func,a,b,ss)
	!use global
	integer, parameter :: jmax = 20, jmaxp = jmax + 1, k = 5, km = k - 1
          
	real (kind = 8) :: a, b, func, ss
	real (kind = 8), parameter :: eps=1.D-6
	external func

!--------------- ---- Usa polint e trapezoidal ---------------------------------
	integer :: j
	real (kind = 8) :: dss, h(jmaxp), s(jmaxp)
	
	h(1) = 1
	do j= 1, jmax
		call trapzd(func, a, b, s(j), j)
		
		if (j >= k) then
			call polint(h(j-km), s(j-km), k, 0.0D+0, ss, dss)
			
			if (abs(dss) <= eps*abs(ss)) then
				exit
			end if
			
		end if
		
		s(j+1) = s(j)
		h(j+1) = 0.25*h(j)
	end do
    
    if (abs(dss) > eps*abs(ss)) then 
		print *, "muitos passos na qromb"
		call exit(1)
	end if
end subroutine qromb1
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine trapzd(func, a, b, s, n)
	!use global
	integer :: n, it, j
	
	real (kind = 8) :: a, b, s, func, del, sum, tnm, x
	external func
	
	if (n == 1) then
		s=0.5*(b-a)*(func(a) + func(b))
	else
		it= 2**(n-2)
		tnm = it
		del = (b-a)/tnm
		x = a+(0.5)*del
		sum = 0.
		do j= 1, it
			sum= sum + func(x)
			x = x + del
		end do
		
		s=0.5*(s + (b-a)*sum/tnm)
	end if
end subroutine trapzd
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine polint(xa, ya, n, x, y, dy)
	!use global
	integer :: n, i, m, ns
	integer, parameter :: nmax = 10
	real (kind = 8) :: dy, x, y, xa(n), ya(n)
	real (kind = 8) :: den, dif, dift, ho, hp, w, c(nmax), d1(nmax)

	ns= 1
	dif = abs(x-xa(1))

	do i = 1, n
		dift = abs(x-xa(i))
		if(dift < dif) then
			ns= i
			dif = dift
		end if
		
		c(i) = ya(i)
		d1(i) = ya(i)
	end do
	
	y = ya(ns)
	ns= ns - 1
	do m=1, n-1
		do i=1, n-m
			ho = xa(i)-x
			hp = xa(i+m) - x
			w = c(i+1) - d1(i)
			den = ho-hp
			if (den == 0.) then
				print*, "falha na polint"
				call exit(1)
			end if
			
			den = w/den
			d1(i)= hp*den
			c(i)= ho*den
		end do
		
		if (2*ns < n-m) then
			dy = c(ns+1)
		else
			dy = d1(ns)
			ns = ns -1
		end if
		
		y = y+dy
	end do
end subroutine polint
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
