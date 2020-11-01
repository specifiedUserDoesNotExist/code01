!-------------------------------------------------------------------------------
!                                                                                 !
!       MODULO DE DISTANCIA SUPERNOVAS MARGINALIZADO, CHI QUADRADO                !
!                                                                                 !
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
module global

    implicit none
    real (kind = 8) :: wm, wl
    
    integer :: NP! = 580
    integer :: SP! = 6
    integer, parameter :: d = 2
    integer :: NDOF! = SP + NP  - d
    integer, parameter :: INT = 100
    integer :: param
    
    real (kind = 8), dimension(:), allocatable :: redz, miobs, sigma
    real (kind = 8), dimension(:), allocatable :: z, Rcb, xsRcb
    
    real (kind = 8), parameter :: PI = 3.1415926536d0
    real (kind = 8) :: A, B
    real (kind = 8), dimension(:), allocatable :: Av
    
end module global
!-------------------------------------------------------------------------------
function prob(c3, n)
    
    use global
    implicit none
    
    integer :: n
    real (kind = 8), intent(in) :: c3
    real (kind = 8) :: prob, gammp
    
    if ( n == 68 ) then
        prob = gammp(NDOF/2.d0, 0.5d0*(c3 + 2.3d0))
    else if ( n == 95 ) then 
        prob = gammp(NDOF/2.d0, 0.5d0*(c3 + 6.17d0))
    else if ( n == 99) then 
        prob = gammp(NDOF/2.d0, 0.5d0*(c3 + 11.8d0))
    else
        prob = gammp(NDOF/2.d0, 0.5d0*c3)
    end if

end function prob
!-------------------------------------------------------------------------------
function func1(x)
    
    use globaL
    real (kind = 8) :: x, func1
    
    !func1 = 1.0d0/sqrt(wm*((1.d0+x)**3.d0)+(1.d0 - wm - wl)*((1.d0+x)**2.d0)+wl)
    func1 = 1.0d0/sqrt(wl + ((1.d0+x)**2.d0)*(1.d0+x*wm-wl))

end function func1
!-------------------------------------------------------------------------------
function Rcbth(wmV, wlV, zk, xint)

	use global
	implicit none
	
	real (kind = 8) :: Rcbth
	real (kind = 8), intent(in) ::  wmV, wlV, zk, xint
	real (kind = 8) :: wk, sqwk, v
	real (kind = 8) :: xxint, nu, de
	
	external func1
	
	wk = 1.d0 - wlV - wmV
	sqwk = sqrt(abs(wk))
	v = 1.d0 + zk
	
	call qromb1(func1, 0.d0, zk, xxint)
	
	if ( wk < 0 ) then
		nu = (sin(xint*sqwk)*((v**2.d0)*(wm*v+wk)+wl)**(1.d0/6.d0))/sqwk
     	de = ((sin(xxint*sqwk))/sqwk)**(2.d0/3.d0)*zk**(1.d0/3.d0)
	else if ( wk == 0 ) then
		nu = xint*((v**2.d0)*(wm*v+wk)+wl)**(1.d0/6.d0) 
		de = (xxint**(2.d0/3.d0))*zk**(1.d0/3.d0)
	else
		nu = (sinh(xint*sqwk)*(((v**2.d0)*(wm*v+wk)+wl))**(1.d0/6.d0))/sqwk 
		de = ((sinh(xxint*sqwk))/sqwk)**(2.d0/3.d0)*(zk)**(1.d0/3.d0)
	end if

	Rcbth = nu/de

end function Rcbth
!-------------------------------------------------------------------------------
function hdL(wmV, wlV, redzk, miobsk)
    
    use global
    implicit none
    
    real (kind = 8) :: hdL
    real (kind = 8), intent(in) :: wmV, wlV, redzk, miobsk
    real (kind = 8) :: wk, sqwk
    real (kind = 8) :: intE, varAux
    
    external func1

    wk = 1.d0 - wlV - wmV
    sqwk = sqrt(abs(wk))
    
    call qromb1(func1, 0.d0, redzk, intE)
    
    if ( wk < 0 ) then 
        varAux = (1+redzk)*sin(sqwk*intE)/sqwk
    else if ( wk == 0 ) then
        varAux = (1+redzk)*intE
    else
        varAux = (1+redzk)*sinh(sqwk*intE)/sqwk
    end if
    
    !alpha
    if ( param == 1 .or. param == 3 .or. param == 5) then
        hdL =  5.d0*log10(varAux) - miobsk
    !magth
    else if ( param == 2 ) then
        hdL =  5.d0*log10(varAux) + 42.384d0 - 5.d0*log10(300.d0) - miobsk
    else if ( param == 4 ) then 
        hdL = 5*log10(varAux) + 42.384d0 - 5*log10(0.738d0) - miobsk    !h=0.738    
    end if
    
end function hdL
!-------------------------------------------------------------------------------
function chi2fun(wmV, wlV)
    
    use global
    implicit none
    
    real (kind = 8) :: C, xint, cmbbao, t1, t2, Rcbth
    real (kind = 8) :: hdlP, varAux, chi2fun, hdL
    real (kind = 8), intent(in) :: wmV, wlV
    
    integer :: k
    
    external func1
    
    B = 0.d0
    C = 0.d0
    
    do k = 1, NP
        hdlP = hdL(wmV, wlV, redz(k), miobs(k))
        varAux = hdlP*Av(k)
        B = B + varAux
        C = C + hdlP*varAux
    end do
    
    if ( param == 1) then
        chi2fun =  C - B*B/A
    else if ( param == 2 ) then 
        chi2fun =  C - (B/A)*(B+0.921034d0) - 2.d0*log(0.678d0)
    else if ( param == 3 ) then 
        chi2fun = C - B*B/A + log(A) - log(2.d0*PI)
    else if ( param == 4 ) then
        chi2fun =  C
	else if ( param == 5 ) then
		
		call qromb1(func1, 0.d0, 1090.d0, xint)
		cmbbao = 0.d0
		
		do k = 1, SP
			t1 = (Rcb(k)-Rcbth(wmV, wlV, z(k), xint))**(2.d0)
			t2 = xsRcb(k)**(2.d0)
		cmbbao = cmbbao + t1/t2
    	end do
		
		chi2fun = C - B*B/A + cmbbao
    end if
    
end function chi2fun
!-------------------------------------------------------------------------------
subroutine invMatrixFisher(wmmin, wlmin, chi2, err1, err2)
    
    !use global
    implicit none
    
    real (kind = 8) :: wmmin, wlmin, chi2fun, dx2, dxdy, dy2, errDx, errDy, errDxDy
    real (kind = 8) :: chi2, err1, err2, varAux!, derivada
    
    call derivada(wmmin, wlmin, chi2, 1, dx2, errDx)
    call derivada(wmmin, wlmin, chi2, 2, dy2, errDy)
    call derivada(wmmin, wlmin, chi2, 3, dxdy, errDxDy)
    
    dx2 = dx2/2.d0
    dy2 = dy2/2.d0
    dxdy = dxdy/2.d0
    
    varAux = 1.d0/(dx2*dy2-dxdy*dxdy)
    err1 = sqrt(dy2*varAux)
    err2 = sqrt(dx2*varAux)

end subroutine invMatrixFisher
!-------------------------------------------------------------------------------
program main

    use global
    implicit none
    
    external func1
    
    integer :: i, m, l
    real (kind = 8) :: SNe, probP, c6min, chi2fun, prob, chi2P
    real (kind = 8) :: wlmin, wmmin, Bmin, xint, cmbbaoF
    real (kind = 8) :: prob68, prob95, prob99, cmbbao
    real (kind = 8), dimension(2) :: error
    
    real (kind = 8) :: c5min, wlminAux, wmminAux
    
    character (len = 30) :: arquivoDeEntradaSN
    character (len = 30) :: arquivoDeEntradaRCB
    character (len = 30) :: arquivoDeSaida
    character (len = 30) :: arquivoSaidaPython
    character (len = 5) :: paramString
    
    integer :: numeroParametro = 1
    integer :: numeroArquivoDeEntradaSN = 2
    integer :: numeroArquivoDeEntradaRCB
    integer :: numeroArquivoDeSaida
    integer :: numeroArquivoSaidaPython
    
    
    call get_command_argument(numeroParametro, paramString)
    read (paramString, *) param
    
    call get_command_argument(numeroArquivoDeEntradaSN, arquivoDeEntradaSN)
    
    if ( param == 5 ) then
    	numeroArquivoDeEntradaRCB = 3
    	numeroArquivoDeSaida = 4
    	numeroArquivoSaidaPython = 5
    	call get_command_argument(numeroArquivoDeEntradaRCB, arquivoDeEntradaRCB)
    else
    	numeroArquivoDeSaida = 3
    	numeroArquivoSaidaPython = 4
    end if
	
    call get_command_argument(numeroArquivoDeSaida, arquivoDeSaida)
    call get_command_argument(numeroArquivoSaidaPython, arquivoSaidaPython)
       
    !open (unit=numeroArquivoDeEntradaSN, file=arquivoDeEntradaSN)
    open (unit=numeroArquivoDeSaida, file=arquivoDeSaida)
    open (unit=numeroArquivoSaidaPython, file=arquivoSaidaPython)    
    !---------------------------------------------------------------------------
    
    !-----------  LER O ARQUIVO DE ENTRADA (A AMOSTRA OBSERVADA) ----------
    call readNumberLines(arquivoDeEntradaSN, NP)
    
    open (unit=numeroArquivoDeEntradaSN, file=arquivoDeEntradaSN)
    allocate(redz(NP), miobs(NP), sigma(NP))
    do i = 1, NP
        read(numeroArquivoDeEntradaSN,*) redz(i), miobs(i), sigma(i)
    end do
    close(numeroArquivoDeEntradaSN)
    
    if (param == 5 ) then
    	call readNumberLines(arquivoDeEntradaRCB, SP)
    
		open (unit=numeroArquivoDeEntradaRCB, file=arquivoDeEntradaRCB)
		allocate(z(SP), Rcb(SP), xsRcb(SP))
		open (unit=numeroArquivoDeEntradaRCB, file=arquivoDeEntradaRCB)
		do i = 1, SP
			read(numeroArquivoDeEntradaRCB,*) z(i), Rcb(i), xsRcb(i)
		end do
		close(numeroArquivoDeEntradaRCB)
	end if
    
    !close(numeroArquivoDeEntradaSN)
    !---------------------------------------------------------------------------
    if ( param /= 5 ) then
    	NDOF = NP - d
    else
    	NDOF = NP + SP - d
    end if
    
    c6min = 10000.00 ! VARIAVEL QUE SERA TESTADA
    
    A = 0.d0
    allocate(Av(NP))
    do i = 1, NP
        Av(i) = 1.d0/sigma(i)**2.d0
        A = A + Av(i)
    end do
    
    do m = 0, INT
        do l = 0, INT
        
            wm = m/(1.d0*INT)
            wl = l/(1.d0*INT)

            chi2P = chi2fun(wm, wl)

            if ( chi2P >= 0 ) then
                probP = prob(chi2P, 0)
                write(numeroArquivoDeSaida,*) wm, wl, probP
                
                if ( chi2P < c6min ) then
                    c6min = chi2P
                    wlmin = wl
                    wmmin = wm
                    Bmin = B
                end if
            else
                write(numeroArquivoDeSaida,*) wm, wl, "NaN"
            end if
        end do
    end do
    
    close(numeroArquivoDeSaida)

    ! - Valor das probabilidades 68.3, 95.4 e 99.7% para desenho curva --
    prob68 = prob(c6min, 68)  !68,3% - 1 sigma
    prob95 = prob(c6min, 95)  !95,4% - 2 sigma
    prob99 = prob(c6min, 99)  !99,7% - 3 sigma

    call invMatrixFisher(wmmin, wlmin, c6min, error(1), error(2))

    write(*,*) param, " -> ", c6min, wmmin, wlmin, c6min
    
    write(numeroArquivoSaidaPython, *) "ndof ", NDOF
    write(numeroArquivoSaidaPython, *) "wm ", wmmin
    write(numeroArquivoSaidaPython, *) "wl ", wlmin
    write(numeroArquivoSaidaPython, *) "wk ", 1.d0 - wlmin - wmmin
    write(numeroArquivoSaidaPython, *) "Mu ", -Bmin/A
    write(numeroArquivoSaidaPython, *) "h ", 10.0d0**((42.384d0+(Bmin/A))/5.0d0)
    write(numeroArquivoSaidaPython, *) "c6min ", c6min
    write(numeroArquivoSaidaPython, *) "chi2 ", c6min/NDOF
    write(numeroArquivoSaidaPython, *) "prob68 ", prob68
    write(numeroArquivoSaidaPython, *) "prob95 ", prob95
    write(numeroArquivoSaidaPython, *) "prob99 ", prob99
    write(numeroArquivoSaidaPython, *) "err1 ", error(1)
    write(numeroArquivoSaidaPython, *) "err2 ", error(2) 
    
    close(numeroArquivoSaidaPython)
  	deallocate(Av, redz, miobs, sigma)
  	if ( param == 5 ) then
  		deallocate(z, Rcb, xsRcb)
  	end if
  
end program main
