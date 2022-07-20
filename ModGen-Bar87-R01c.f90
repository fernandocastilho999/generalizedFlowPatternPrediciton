





INTEGER FUNCTION ModelGBar87(D, g, DenL, DenG, VisL, VisG, jL, jG, sTen, AngR)


    IMPLICIT NONE
    DOUBLE PRECISION :: jL, jG, AngR
    DOUBLE PRECISION :: D, g, DenL, DenG, VisL, VisG, sTen
    DOUBLE PRECISION :: Crit1a, Crit2a,Crit3a, Crit4a, Crit4b, &
            Crit5a, Crit5b, Crit6a, Crit6b, Crit6c, Crit7a

    DOUBLE PRECISION :: ReL, ReG, X2, X, Y, F, s, K, T, Z,W
    DOUBLE PRECISION :: A_L, A_G, A_P, u_L, u_G, S_I, S_L, D_L
    DOUBLE PRECISION :: CL, CG, n, m
    DOUBLE PRECISION :: fM, dC, dCD, dCB
    DOUBLE PRECISION :: alfaL, Calc_alfaL, Rsm
    DOUBLE PRECISION :: Uo, gam, Clift
    DOUBLE PRECISION :: alfaC, alfaS, Rs

    DOUBLE PRECISION :: h_L, Calc_h_L2, h_x


    DOUBLE PRECISION, PARAMETER :: PI = 3.145159265



    !PRINT*,'QUESTION1: There is Bubble Flow?'


    ModelGBar87 = 0;

    ! Reynolds numbers
    ReL = DenL * jL * D / VisL;
    ReG = DenG * jG * D / VisG;
    ! Determination of CL - CG and n-m exponents of PressureDrop equation
             !IF (ReL<4100) THEN
                CL=16.
                n = 1.
            !ELSE
                CL = 0.046
                n = 0.2
            !ENDIF

!            IF (ReG<4100) THEN
!                CG=16.
!                m = 1.
!            ELSE
                CG = 0.046
                m = 0.2
           ! ENDIF

    ! Determination of
    !!!!!!!!

        fM = CL * ( (jL+jG) * D / (VisL/DenL) )**-n


        dCD = 2. * ( 0.4*sTen / ( (DenL-DenG)* g ) )**0.5

        dCB = 3./8. * (DenL / (DenL-DenG)) * (fM * (jL+jG)**2.) / (g * MAX(1e-3,COS(AngR)))

        dC = MIN(dCB,dCD)

        Crit1a = ( 0.725 + 4.15 * (jG/(jL+jG))**0.5 ) * ( sTen/DenL )**(3./5.) * ( 2*fM/D * (jL+jG)**3 )**(-2./5.)

        !Print*,'Crit1',Crit1

   !!!!!!!!
    ! Parameter of Lockhart-Martinelli
    X2 = ( (4.*CL/D) * ( DenL*jL*D / VisL )**-n * (DenL*(jL**2.)/2.) ) / &
        ( (4.*CG/D) * ( DenG*jG*D / VisG )**-m * (DenG*(jG**2.)/2.) )

    !print*, X2
    X = SQRT(X2)

    Y = ( (DenL-DenG)*g*SIN(AngR) ) / ( (4.*CG/D) * ( DenG*jG*D / VisG )**-m * (DenG*(jG**2)/2) )


    ! Obtain hL
    h_L=Calc_h_L2(X2, Y, m , n, jL, jG)
    !PRINT*, 'jL=',jL,'jG=',jG,'X=',X,'Y=',Y,'h_L=',h_L


    ! Calculate F
    F = SQRT( (DenG)/(DenL-DenG) ) * (jG) / SQRT(D * g * COS(AngR))

    !Criterion*

    h_x = 2. * h_L - 1.
    A_L = 0.25 * ( PI - ACOS(h_x) + h_x * SQRT(1 - h_x**2.) )
    A_G = 0.25 * ( ACOS(h_x) - h_x * SQRT(1 - h_x**2.) )
    A_P = A_L + A_G
    u_L = A_P / A_L
    u_G = A_P / A_G

    S_I = SQRT( 1 - h_x**2 )

    S_L = PI - ACOS(h_x)
    D_L = 4 * A_L / S_L

    Crit2a = F**2. * ( 1./((1-h_L)**.2) * (u_G**2. * S_I)/(A_G))


    Z = ( (4.*CL/D) * ( DenL*jL*D / VisL )**(-n) * (DenL*(jL**2)/2.) ) /&
            ( DenL * g * COS(AngR) )
    Crit3a = 2.* (A_L/A_P)**2. * (1.-h_L)


    s = 0.01
    K = F * SQRT(ReL)
    Crit4a = 2. / ( (SQRT(u_L))*(u_G)*(SQRT(s)) );

    W = jL / SQRT(g*D)
    Crit4b = 1.5 * SQRT(h_L) * ( A_L / A_P )


    ! ANULAR OR NOT!
    alfaL = Calc_alfaL(X2,Y)
    !print*,jL,jG,'alfaL = ', alfaL
    Crit5a = ( (2. - 3./2.*alfaL) / ( (alfaL**3.) * (1.-3./2.*alfaL) ) ) * X2

    Rsm = 0.48
    Crit5b = alfaL/Rsm

    ! THERE IS BUBBLE?
    Crit6a = 19. * ( ( ( (DenL-DenG)*sTen ) / ( (DenL**2)*g ) )**(1./2.) )

    Clift = 0.8
    gam = 1.3
    Uo = 1.53 * ( ( ( g*(DenL-DenG)*sTen ) / (DenL**2.) )**(1./4.) )
    Crit6b = (3./4.) * COS(PI/4.) * ( (Uo**2.)/g ) * ( Clift*(gam**2.)/D )
    Crit6c = COS(AngR) / (SIN(AngR))**2.

    ! BUBBLE OR INTERMITTENT
    alfaC = 0.25

    Crit7a = ((1.-alfaC)/alfaC)*jG - 1.53 * (1.-alfaC) * &
        ( ( ( g*(DenL-DenG)*sTen ) / (DenL**2.) )**(1./4.) ) * SIN(AngR)


    ! SLUG OR ELONGATED!
    alfaS = 0.058 * ( ( dC * ( ( (2.*fM/D) * ((jL+jG)**3.) )**(2./5.) ) * ( (DenL/sTen)**(3./5.) ) - 0.725 )**2. )
    Rs = 1. - alfaS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHOOSE MODELS!!
    IF ( (dC.GE.Crit1a) .AND. (jL.GE.(jG*0.48/0.52)) ) THEN !!!!!FECHADO1!!!!

        ModelGBar87 = 20

    ELSEIF ( (Crit2a.LT.1.) .AND. ( Z.LT.Crit3a ) .AND. &
         (( (AngR.LT.0.) .AND. (W.LT.Crit4b) ) .OR. &
        ( (AngR.GE.0.) .AND. (K.LT.Crit4a) )) ) THEN

        ModelGBar87 = 21

    ELSEIF ( (Crit2a.LT.1.) .AND. ( Z.LT.Crit3a ) ) THEN

        ModelGBar87 = 22

    ELSEIF ( (Y.LT.Crit5a) .AND. (Crit5b.LT.0.5) ) THEN

        ModelGBar87 = 23

    !(AngR.GE.(80.*PI/180.))
    ELSEIF ( ( (D.LE.Crit6a) ) .AND. &
            ( ( (jL.GT.Crit7a) .AND. (AngR.GE.(60.*PI/180.)) ) ) ) THEN

        ModelGBar87 = 24

    ELSEIF ( (Rs.LT.1) .AND. (AngR.EQ.0.) .AND. ( Rs.GT.0.48 ) ) THEN

        ModelGBar87 = 27

    !ELSEIF ( Rs.LE.0.48 ) THEN

     !   ModelGBar87 = 25

    ELSEIF ( (Rs.LT.0.48).AND.(AngR.GT.0.) .AND. (AngR.GE.(70.*PI/180.))  ) THEN

        ModelGBar87 = 26




    ELSE
        ModelGBar87 = 25
    ENDIF
    RETURN


END


DOUBLE PRECISION FUNCTION Calc_h_L2(X2, Y, m , n, jL, jG)
    IMPLICIT NONE
    DOUBLE PRECISION :: jL, jG
    DOUBLE PRECISION ::  f, f1, f2
    DOUBLE PRECISION :: m, n
    DOUBLE PRECISION :: h_L, h_L1, h_L2
    DOUBLE PRECISION:: X2, Y
 !   INTEGER :: nx
    DOUBLE PRECISION :: fMSec2
    DOUBLE PRECISION :: tol, var

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    h_L1 = 1e-5
    h_L2 = 1. - 1.e-5

    tol = 1e-7
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        var = 1.
        DO WHILE (var.GE.tol)
            f1 = fMSec2(h_L1, X2, Y, m, n)
            f2 = fMSec2(h_L2, X2, Y, m, n)

            h_L = (h_L1 + h_L2) / 2.
            f = fMSec2(h_L, X2, Y, m, n)

            IF (f*f1 .LT. 0.) THEN
                h_L2 = h_L
            ELSE
                h_L1 = h_L
            ENDIF

            var = (ABS(h_L1-h_L2))

        ENDDO

    Calc_h_L2 = h_L
    RETURN

END


DOUBLE PRECISION FUNCTION fMSec2(h_L, X2, Y, m, n)
    IMPLICIT NONE
    DOUBLE PRECISION :: h_L, h_x
    DOUBLE PRECISION :: m, n, X2, Y
    DOUBLE PRECISION :: A_L, A_G, S_L, S_G, S_I, A_P, u_L, u_G, D_L, D_G
    DOUBLE PRECISION, PARAMETER :: PI = 3.145159265


    h_x = 2.*h_L-1.


    A_L = 0.25 * ( PI - ACOS(h_x) + h_x * SQRT(1. - h_x**2.) )
    A_G = 0.25 * ( ACOS(h_x) - h_x * SQRT(1. - h_x**2.) )
    S_L = PI - ACOS(h_x)
    S_G = ACOS(h_x)
    S_I = SQRT( 1. - h_x**2. )
    A_P = A_L + A_G
    u_L = A_P / A_L
    u_G = A_P / A_G
    D_L = 4. * A_L / S_L
    D_G = 4. * A_G / (S_G + S_I)

    fMSec2 = (X2**1.) * ( ((u_L*D_L)**(-n)) * (u_L**2.) * (S_L / A_L) ) -&
        ( ((u_G*D_G)**(-m)) * (u_G**2.) * (S_G/A_G + S_I/A_L + S_I/A_G) ) + 4.*Y

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION Calc_alfaL(X2,Y)

IMPLICIT NONE
DOUBLE PRECISION :: fMBis2
DOUBLE PRECISION :: X2, Y, alfaL1, alfaL2, alfaL, f1, f2, f
DOUBLE PRECISION :: tol, var

alfaL1 = 1e-4
alfaL2 = 1. - 1e-4

tol = 1e-7
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        var = 1.
        DO WHILE (var.GE.tol)
            f1 = fMBis2(alfaL1,X2,Y)
            f2 = fMBis2(alfaL2,X2,Y)

            alfaL = (alfaL1 + alfaL2) / 2.
            f = fMBis2(alfaL,X2,Y)

            IF (f*f1 .LT. 0.) THEN
                alfaL2 = alfaL
            ELSE
                alfaL1 = alfaL
            ENDIF

            var = (ABS(alfaL1-alfaL2))

        ENDDO

    Calc_alfaL = alfaL
    RETURN

END


DOUBLE PRECISION FUNCTION fMBis2(alfaL,X2,Y)

IMPLICIT NONE
DOUBLE PRECISION :: alfaL, X2, Y


    fMBis2 = ( (1.+75.*alfaL) / ( ((1.-alfaL)**(5./2.))*alfaL ) ) - ( (1./(alfaL**3.))*X2 ) - Y


END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

