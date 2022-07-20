! #################################################################
! NFPC R0.1
! Program "NUEM Flow Pattern Calculator"

! ##################################################################
! This program has been developed by Fernando Castillo Vicencio
! at the NUEM-UTFPR, with the aim to predict the Flow Pattern that
! occurs for every JL and JG.
! ##################################################################
! mainR01a is the main program of the project
PROGRAM mainR01a

    ! #################################################
! (1) Declaration of Variables
    IMPLICIT NONE

    ! ------------------------------
    ! (1.1) Constants
        DOUBLE PRECISION  ::  D, Ang, g, Po, ToC, L
            ! D : Diameter of Pipe (m)
            ! g : Acceleration of Gravity (m/s²)
            ! Po : Ambient Pressure (Pa)
            ! ToC : Ambient Temperature (°C)

    ! ------------------------------
    ! (1.2) ExecutionVariables
        DOUBLE PRECISION        ::  DenL, DenG, VisL, VisG, sTen
        INTEGER     ::  ModelTai80, ModelHTai76, ModelGBar87
        DOUBLE PRECISION :: jL, jG, AngR
            ! jL   : Superficial velocity of Water (m/s)
            ! jG   : Superficial velocity of Air (m/s)
            ! DenL : Density of Water (kg/m³)
            ! DenG : Density of Air (kg/m³)
            ! VisL : Viscosity of Water (Pa.s)
            ! VisG : Viscosity of Air (Pa.s)
            ! sTen : Surface Tension of Water (N/m)


    ! ------------------------------
    ! (1.3) AuxiliarVariables
        DOUBLE PRECISION  :: a0, an, To
        DOUBLE PRECISION :: fPat
        DOUBLE PRECISION, PARAMETER :: PI = 3.145159265
        INTEGER :: i, j, n

            ! Crit : Criterion to evaluate flow patterns
            ! To : Ambient Temperature (K)

    ! ##################################################
    ! (2) Data

    ! ------------------------------
    !   (2.1) Operation Data
        !jL = 2
        !jG = 0.5
        D = 0.051
        Ang=80.
        L = 9.

        g = 9.81
        Po = 101325.
        ToC = 25.

    ! ------------------------------
    !   (2.2) Auxiliar Variables
        To = ToC + 273.15
        AngR = Ang * PI / 180

    ! ##################################################
    ! (3) Execution

    ! ------------------------------
    !   (3.1) Properties of Fluids

        ! Density of Water as function of Temperature (°C) by
        ! Tanaka, M., et al. "Recommended table for the density of water between 0 C and 40 C based on recent experimental reports." Metrologia 38.4 (2001): 301.
        DenL = 999.974950 * ( 1 - (ToC-3.983035)**2 * (ToC+301.797) / (522528.9*(ToC+69.34881)) )
        PRINT*, "WaterDensity:", DenL, "kg/m³"



        ! Viscosity of Water as function of Temperature (Pa.s) by
        ! Kestin, Joseph, Mordechai Sokolov, and William A. Wakeham. "Viscosity of liquid water in the range− 8 C to 150 C." Journal of Physical and Chemical Reference Data 7.3 (1978): 941-948
        VisL = 1002.e-6 * EXP ( (20.-ToC)/(ToC+95) * ( 1.2378 - 1.303e-3*(20.-ToC) + 3.06e-6*(20.-ToC)**2 + 2.55e-8*(20.-ToC)**3 ) )
        PRINT*, "WaterViscosity:", VisL, "Pa.s"

       ! Density of Air by Gas Ideal law
        DenG = Po / (287. * To)
        PRINT*, "AirDensity:", DenG, "kg/m³"

        ! Viscosity of Air as function of Temperature (Pa.s) by
        ! Sutherland, William. "LII. The viscosity of gases and molecular force." The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science 36.223 (1893): 507-531.
        VisG = 1.512041288e-6 * To**1.5 / (To + 120.)
        PRINT*, "AirViscosity:", VisG, "Pa.s"

        ! Surface Tension of Water as function of Temperature (N/m) by
        ! Vargaftik, N. B., B. N. Volkov, and L. D. Voljak. "International tables of the surface tension of water." Journal of Physical and Chemical Reference Data 12.3 (1983): 817-820.
        sTen = 235.8e-3 * ( (647.15-To) / (647.15) )**1.1256 * ( 1 - 0.625 * (647.15-To) / (647.15) )
        PRINT*, "SurfaceTension:", sTen, "N/m"


    ! ------------------------------
    !   (3.2) There is Bubble Flow?
       ! Crit1 = ( (DenL**2)*g*(D**2) ) / ( (DenL-DenG)*(sTen**0.25) )
       ! PRINT*, 'Crit1=',Crit1
        !fPat = ModelTai80(D, g, Po, ToC, L, DenL, DenG, VisL, VisG, sTen, jL, jG)

        OPEN (10, file='Output_R01a_testMarco.txt', status='unknown')
        n=30.

        DO i=1,n
            a0=0.01 ; an=50;
            jL = a0 * (an/a0)**((i-1.)/(n-1.))


            DO j=1,n
                a0=0.01; an=80;
                jG = a0 * (an/a0)**((j-1.)/(n-1.))


                !IF (Ang.EQ.0.) THEN
                    !PRINT*,'HorizontalFlow'
                    !jL = 0.1; jG = 3.;
                    !D=0.05;
                    !DenL=993.; VisL = 6.8e-4; DenG = 1.14; VisG=1.9e-5
                    !fPat = ModelHTai76(D, g, DenL, DenG, VisL, VisG, jL, jG, AngR)

                !ELSEIF (Ang.EQ.90) THEN

                    !PRINT*,'VerticalFlow'


                 !   fPat = ModelGBar87(D, g, L, DenL, DenG, VisL, sTen, jL, jG, AngR)

                !ELSEIF ( (Ang.LT.90.) .AND. (Ang.GT.-90.) ) THEN
                  !PRINT*,'InclinatedFlow'
                   fPat = ModelGBar87(D, g, DenL, DenG, VisL, VisG, jL, jG, sTen, AngR)


                !ELSE

                !    PRINT*,'TA ERRADO!!'

                !ENDIF



                WRITE(10,*) jL,   jG, fPat
            ENDDO


        ENDDO

        CLOSE(10)


   end program mainR01a
