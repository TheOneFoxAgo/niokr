module constants
    Real(8) :: Pi=acos(-1.0d0)
    Real(8) :: TC, Gamma, gammaspr, gammar, c, h, T, AMASS, ATEMASS, omegahfs, omegar1, omegar2, Deltaomegar, Massa, &
               Vt, Vtr, Omega1, Omega2, V10R, V20R, V10I, V20I, k1, k2, Omega1P, Omega2P, Omega1M, Omega2M, WP,     &
               DeltaP, WM, DeltaM, GammaL, d31, d32, na, a
    Real(8) :: Vzleft, Vzright
    Real(8) :: alpha, betta, Nlev, factor
  contains
   subroutine get_consants()
    implicit none

    ! PARAMETRY YACHEIKI i LASERA
    TC          = 50.0d0                     ! (* temperatura, Celsium *)
    Gamma       = 100.0d0          
    gammaspr    = 2.0d0*Pi*6.0d6          ! (* postoyannaya spontannogo raspada, s^-1 *)
    gammar      = 1.06d0*gammaspr  !  
    
    ! CONSTANTY I PREOBRASOVANIE PARAMETROV K EDINITSAM, ISPOLZUEMYM V RASCHETE
    
!                                         (* relaxaciya dla opticheskix kogerentciy, s^-1 *)
    c          = 2.997d10;   !  2.9979245800d10        ! (* scorost sveta, cm/s *)
    h          = 1.0546d-27;                           ! (* postoyannaya Planka, erg*s *)
    T          = (TC+273.0d0)/11604.5d0;               ! (* temperatura, eV *)
    AMASS      = 86.90918d0;                       ! (* atomnaya massa *)
    ATEMASS    = 9.315d8;                        ! (* edinica atomnoy massy, eV *)
    
    omegahfs   = 2.0d0*Pi*6834.0d0*1.0d6
!    omegahfs   = 0.0d0
    omegar1    = 2.372d15;                         ! (* chastota perexoda D2, s^-1 *)
    Deltaomegar = omegahfs
    omegar2    = omegar1 - Deltaomegar
    Massa      = AMASS*ATEMASS                     ! (* massa atoma, eV *)
    Vt         = Sqrt((2.0d0*T)/Massa)             ! (* teplovaya skorost *)
    VTr        = Vt*c                              ! (*otnositelnaya skorost atomov,cm/s *)
    GammaL     = Gamma*2.0d0 ! Gamma/2.0d0 !  
    k1         = omegar1/c 
    k2         = omegar2/c  
!    d31        = 0.0d0 !
    d31        = sqrt(3.0d0*gammaspr*h*c**3/(2.0d0*omegar1**3))
!    d32        = 0.0d0 !
    d32        = sqrt(3.0d0*gammaspr*h*c**3/(2.0d0*omegar2**3))
!    na         = 4.0d10 ! 7.0d10 ! 2.5d10 ! 1.0d10 ! 0.8d11 ! 0.4d11 ! 0.2d11 ! 1.31d11 ! 1.0d11 ! 3.0d11
    na         = exp(-3.0d-4*Tc**2+0.119d0*Tc+20.4d0)
    Vzleft     = 0.0d0     ! 10.0d0    ! 
    Vzright    = 6.0d0*Vtr ! 3.0d0*Vtr !

    Omega1     = 0.0d0  
    
    Omega2     = 0.0d0
    V10R       = 5.0d5    
    V20R       = 5.0d5
    
    V10I       = 0.0d0 ! 1.0d4 ! 
    V20I       = 0.0d0 ! 1.0d4 !
    
    a          = 1.1d0 ! 2.20d0 
    
    alpha      = 0.0d0 ! 1.0d0   
    betta      = 0.0d0 ! 1.0d0  
    Nlev       = 2.0d0
    
   end subroutine get_consants

   subroutine write_consants()
    implicit none
    open(1,file='parameters.dat')
    ! PARAMETRY YACHEIKI i LASERA
    write(1,"(' TC         =',es14.6)") TC                     ! (* temperatura, Celsium *)
    write(1,"(' Gamma      =',es14.6)") Gamma
    write(1,"(' gammaspr   =',es14.6)") gammaspr               ! (* postoyannaya spontannogo raspada, s^-1 *)
    write(1,"(' gammar     =',es14.6)") gammar
    
    ! CONSTANTY I PREOBRASOVANIE PARAMETROV K EDINITSAM, ISPOLZUEMYM V RASCHETE
    
!                                         (* relaxaciya dla opticheskix kogerentciy, s^-1 *)
    write(1,"(' c          =',es14.6)") c                              ! (* scorost sveta, cm/s *)
    write(1,"(' h          =',es14.6)") h                              ! (* postoyannaya Planka, erg*s *)
    write(1,"(' AMASS      =',es14.6)") AMASS                          ! (* atomnaya massa *)
    write(1,"(' ATEMASS    =',es14.6)") ATEMASS                        ! (* edinica atomnoy massy, eV *)
    write(1,"(' omegahfs   =',es14.6)") omegahfs
    write(1,"(' omegar1    =',es14.6)") omegar1                        ! (* chastota perexoda D2, s^-1 *)
    write(1,"(' Deltaomegar=',es14.6)") Deltaomegar                    ! (* raznost' chastot perexodov, s^-1 *)
    write(1,"(' omegar2    =',es14.6)") omegar2
    write(1,"(' Massa      =',es14.6)") Massa                          ! (* massa atoma, eV *)
    write(1,"(' Vt         =',es14.6)") Vt                             ! (* teplovaya skorost *)
    write(1,"(' VTr        =',es14.6)") Vtr                            ! (*otnositelnaya skorost atomov,cm/s *)
    write(1,"(' GammaL     =',es14.6)") GammaL
    write(1,"(' k1         =',es14.6)") k1
    write(1,"(' k2         =',es14.6)") k2
    write(1,"(' d31        =',es14.6)") d31
    write(1,"(' d32        =',es14.6)") d32
    write(1,"(' na         =',es14.6)") na
    write(1,"(' T, [eV]    =',es14.6)") T                              ! (* temperatura, eV *)
    write(1,"(' 2.0*Pi*na*k1*d31**2/h =',es14.6)") 2.0d0*Pi*na*k1*d31**2/h
    write(1,*) 
   
    write(1,*) 
    
    write(1,"(' Omega1     =',es14.6)") Omega1
    write(1,"(' Omega2     =',es14.6)") Omega2
    write(1,"(' Re(V1(0))  =',es14.6)") V10R
    write(1,"(' Re(V2(0))  =',es14.6)") V20R
    write(1,"(' Im(V1(0))  =',es14.6)") V10I
    write(1,"(' Im(V2(0))  =',es14.6)") V20I
    write(1,"(' a          =',es14.6)") a
    write(1,"(' Vzleft     =',es14.6)") Vzleft
    write(1,"(' Vzright    =',es14.6)") Vzright
    write(1,*) 
    write(1,"(' alpha      =',es14.6)") alpha
    write(1,"(' betta      =',es14.6)") betta
    write(1,"(' Nlev       =',es14.6)") Nlev

    close(1)
   end subroutine write_consants

end module constants
