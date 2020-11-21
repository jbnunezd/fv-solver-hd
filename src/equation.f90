!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE RiemannSolver
  MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE EvaluateFlux1D
  MODULE PROCEDURE EvaluateFlux1D
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: KappaP1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                :: Prim(1:nVar)
REAL                :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL                :: delta_rho, delta_u, delta_v
REAL                :: x0, y0, xt, xm, ym
REAL                :: sigma, a_shear, v_shear, amplitude, rho0, rho1
REAL                :: rho_postshock, p_postshock, vx_postshock, vy_postshock
REAL                :: rho_preshock, p_preshock, vx_preshock, vy_preshock
REAL                :: MachNumber_shock, MachNumber_vortex
REAL                :: xc(2), xs, a, b, r, d
REAL                :: vm, vtheta, sin_theta, cos_theta
REAL                :: Temp, Temp0, Ta, Tb, Tdiff, Tdiff2
CHARACTER(LEN=255)  :: ErrorMessage
!-------------------------------------------------------------------------------!

Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
  !----------------------------------------------------------------------!
  ! [200] Constant State                                                 !
  !----------------------------------------------------------------------!
  CASE(200)
    Prim(1:nVar) = PrimRefState1(1:nVar)

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [211] Riemann Problem                                                !
  !----------------------------------------------------------------------!
  CASE(211)
    xm = MESH_X0(1)+0.5*MESH_SX(1)
    ym = MESH_X0(2)+0.5*MESH_SX(2)
    IF      ((x(1) .GE. xm) .AND. (x(2) .GE. ym)) THEN
      Prim(1:nVar) = PrimRefState1(1:nVar)
    ELSE IF ((x(1) .LT. xm) .AND. (x(2) .GE. ym)) THEN
      Prim(1:nVar) = PrimRefState2(1:nVar)
    ELSE IF ((x(1) .LT. xm) .AND. (x(2) .LT. ym)) THEN
      Prim(1:nVar) = PrimRefState3(1:nVar)
    ELSE IF ((x(1) .GE. xm) .AND. (x(2) .LT. ym)) THEN
      Prim(1:nVar) = PrimRefState4(1:nVar)
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [214] Double Mach Reflection Problem                                 !
  !----------------------------------------------------------------------!
  CASE(214)
    Prim_in(1:nVar)  = PrimRefState1(1:nVar)
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    x0 = 1.0/6.0
    xt = x0 + (1.0/SQRT(3.0))*x(2)
    IF (x(1) .LT. xt) THEN
      Prim(1:nVar) = Prim_in(1:nVar)
    ELSE
      Prim(1:nVar) = Prim_out(1:nVar)
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [215] Implosion Problem                                              !
  !----------------------------------------------------------------------!
  CASE(215)
    d = 0.15
    IF ((x(1)+x(2)) .GT. d) THEN
      Prim(1:nVar) = PrimRefState1(1:nVar)
    ELSE
      Prim(1:nVar) = PrimRefState2(1:nVar)
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [216] Shock-Vortex Interaction Problem                               !
  !----------------------------------------------------------------------!
  CASE(216)
    xc(1) = 0.25
    xc(2) = 0.50
    xs    = 0.50
    a     = 0.075
    b     = 0.175
    
    MachNumber_shock  = 7.0
    MachNumber_vortex = 1.0
    
    rho_postshock = 1.0
    p_postshock   = 1.0
    vx_postshock  = MachNumber_shock*SQRT(Kappa*p_postshock/rho_postshock)
    vy_postshock  = 0.0
    
    rho_preshock  = rho_postshock*(KappaP1*MachNumber_shock**2) &
                  /(KappaM1*MachNumber_shock**2+2.0)
    p_preshock    = p_postshock*(2.0*Kappa*MachNumber_shock**2-KappaM1)/KappaP1
    vx_preshock   = vx_postshock*rho_postshock/rho_preshock
    vy_preshock   = 0.0
    
    r = SQRT((x(1)-xc(1))**2+(x(2)-xc(2))**2)
    IF (r .LE. b) THEN
      Temp0 = p_postshock/rho_postshock
      vm = MachNumber_vortex*SQRT(Kappa*p_postshock/rho_postshock)
      Ta = (KappaM1*(a**2)*(vm**2))/(Kappa*a**2)/2.0
      Tb = (KappaM1*(a**2))*(-2.0*(b**2)*LOG(a) &
         + (a**2)/2.0-(b**4)/(a**2)/2.0)*(vm**2)/(Kappa*(a**2-b**2)**2)
      Tdiff  = Tb-Ta
      Tdiff2 = (a**2)*KappaM1*(-2.0*(b**2)*LOG(b) &
             + (b**2)/2.0-(b**4)/(b**2)/2.0)*(vm**2)/(Kappa*(a**2-b**2)**2)
      Tdiff2 = Temp0-Tdiff2
      IF (r .LE. a) THEN
        vtheta = vm*r/a
        Temp = KappaM1*(r**2)*(vm**2)/(Kappa*(a**2))/2.0 + Tdiff + Tdiff2
      ELSE IF (r.LE.b) THEN
        vtheta = vm*a/(a**2-b**2)*(r-b**2/r)
        Temp   = KappaM1*(a**2)*(-2.0*(b**2)*LOG(r) + (r**2)/2.0 &
               - (b**4)/(r**2)/2.0)*(vm**2)/(Kappa*(a**2-b**2)**2)+Tdiff2
      END IF
      sin_theta = (x(2)-xc(2))/r
      cos_theta = (x(1)-xc(1))/r
      
      rho_postshock = rho_postshock*(Temp/Temp0)**(Kappa*sKappaM1) 
      p_postshock   = p_postshock*(Temp/Temp0)**(sKappaM1)
      vx_postshock  = vx_postshock - vtheta*sin_theta
      vy_postshock  = vy_postshock + vtheta*cos_theta
    END IF

    IF (x(1) .GT. xs) THEN
      Prim(1) = rho_preshock
      Prim(2) = vx_preshock
      Prim(3) = vy_preshock
      Prim(4) = p_preshock
    ELSE
      Prim(1) = rho_postshock
      Prim(2) = vx_postshock
      Prim(3) = vy_postshock
      Prim(4) = p_postshock
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [217] Kelvin-Helmholtz Instability                                   !
  !----------------------------------------------------------------------!
  CASE(217)
    sigma     = 0.10
    a_shear   = 0.01
    v_shear   = 0.50
    amplitude = 0.10
    rho0      = 0.505
    rho1      = 0.495

    IF (x(2) .GT. 0.0) THEN
      delta_rho = rho0 + rho1*TANH((x(2)-0.5)/a_shear)
    ELSEIF (x(2) .LE. 0.0) THEN
      delta_rho = rho0 - rho1*TANH((x(2)+0.5)/a_shear)
    END IF
    IF (x(2) .GT. 0.0) THEN
      delta_u = +v_shear*TANH((x(2)-0.5)/a_shear)
    ELSEIF (x(2) .LE. 0.0) THEN
      delta_u = -v_shear*TANH((x(2)+0.5)/a_shear)
    END IF
    IF (x(2) .GT. 0.0) THEN
      delta_v = +amplitude*v_shear*SIN(2.0*PI*x(1))*EXP(-(x(2)-0.5)**2/sigma)
    ELSEIF (x(2) .LE. 0.0) THEN
      delta_v = -amplitude*v_shear*SIN(2.0*PI*x(1))*EXP(-(x(2)+0.5)**2/sigma)
    END IF

    Prim(1) = delta_rho
    Prim(2) = delta_u
    Prim(3) = delta_v
    Prim(4) = 1.0

    CALL PrimToCons(Prim,Cons)
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj
!-------------------------------------------------------------------------------!

S = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    ErrorMessage = "Boundary condition defined only for faces 1 and 3"
    WRITE(*,*) ErrorMessage
    STOP
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    ErrorMessage = "Boundary condition defined only for faces 1 and 3"
    WRITE(*,*) ErrorMessage
    STOP
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    Prim_in(1:nVar)  = PrimRefState1(1:nVar)
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        x0 = 1.0/6.0
        xc = MeshBary(1,ii,nElemsY)
        xt = x0 + (1.0+2.0*10.0*t)/SQRT(3.0)
        IF (xc .LT. xt) THEN
          U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
        ELSE
          U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(11) ! Double Mach Reflection
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        x0 = 1.0/6.0
        xc = MeshBary(1,ii,1)
        IF (xc .LT. x0) THEN
          U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
        ELSE
          U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
          U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
        END IF
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_TIMESTEP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds1D(Prim,slowest,fastest)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)           :: Prim(1:nVar)
REAL,INTENT(OUT),OPTIONAL :: slowest
REAL,INTENT(OUT),OPTIONAL :: fastest
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                      :: rho, vx, vy, p, c
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

c  = EOS_SoundSpeed(rho,p)

IF(PRESENT(slowest)) THEN
  slowest = vx-c
END IF

IF(PRESENT(fastest)) THEN
  fastest = vx+c
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: fastestx
REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, p, c
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

c = EOS_SoundSpeed(rho,p)

fastestx = vx+c
fastesty = vy+c

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION EOS_SoundSpeed(rho,p)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: rho, p
REAL            :: EOS_SoundSpeed
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

EOS_SoundSpeed = ABS(Kappa*p/rho)
EOS_SoundSpeed = SQRT(EOS_SoundSpeed)

!-------------------------------------------------------------------------------!
END FUNCTION EOS_SoundSpeed
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_ENERGY, MIN_MOMENTUM
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, Sx, Sy, E, S2
!-------------------------------------------------------------------------------!

rho = Cons(1)
Sx  = Cons(2)
Sy  = Cons(3)
E   = Cons(4)

S2 = Sx*Sx + Sy*Sy

IF (rho .LT. MIN_DENSITY) THEN
  rho = MIN_DENSITY
END IF
IF (E .LT. MIN_ENERGY) THEN
  E = MIN_ENERGY
END IF
IF (ABS(Sx) .LT. MIN_MOMENTUM) THEN
  Sx = 0.0
END IF
IF (ABS(Sy) .LT. MIN_MOMENTUM) THEN
  Sy = 0.0
END IF

Prim(1) = rho
Prim(2) = Sx/rho
Prim(3) = Sy/rho
Prim(4) = (Kappa-1.0)*(E-0.5*S2/rho)


!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_PRESSURE, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, p, v2
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

v2 = vx*vx + vy*vy

IF (rho .LT. MIN_DENSITY) THEN
  rho = MIN_DENSITY
END IF
IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
END IF
IF (ABS(vx) .LT. MIN_SPEED) THEN
  vx = 0.0
END IF
IF (ABS(vy) .LT. MIN_SPEED) THEN
  vy = 0.0
END IF

Cons(1) = rho
Cons(2) = rho*vx
Cons(3) = rho*vy
Cons(4) = (1.0/(Kappa-1.0))*p + 0.5*rho*v2

!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa, sKappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, p, v2
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
p   = Prim(4)

v2  = vx*vx + vy*vy

Flux(1) = rho*vx
Flux(2) = rho*vx*vx + p
Flux(3) = rho*vx*vy
Flux(4) = (Kappa*sKappaM1*p + 0.5*rho*v2)*vx

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolver(PrimL,PrimR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: NormVect(1:nDims,1:nGPs)
REAL,INTENT(IN)  :: TangVect(1:nDims,1:nGPs)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
REAL             :: ConsLL(1:nVar,1:nGPs), ConsRR(1:nVar,1:nGPs)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1,iGP)*PrimL(2,iGP) + NormVect(2,iGP)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1,iGP)*PrimL(2,iGP) + TangVect(2,iGP)*PrimL(3,iGP)
  PrimLL(4,iGP) = PrimL(4,iGP)

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1,iGP)*PrimR(2,iGP) + NormVect(2,iGP)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1,iGP)*PrimR(2,iGP) + TangVect(2,iGP)*PrimR(3,iGP)
  PrimRR(4,iGP) = PrimR(4,iGP)

  CALL PrimToCons(PrimLL(1:nVar,iGP),ConsLL(1:nVar,iGP))
  CALL PrimToCons(PrimRR(1:nVar,iGP),ConsRR(1:nVar,iGP))  

  CALL RiemannSolverByRusanov(&
    ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
    PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))

  ! Rotating back the vector quantities  !
  Flux(2:3,iGP) = NormVect(1:nDims,iGP)*Flux(2,iGP) &
                + TangVect(1:nDims,iGP)*Flux(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolver
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: LambdaMax, fastestL, fastestR
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)
CALL WaveSpeeds1D(PrimL,fastest=fastestL)
CALL WaveSpeeds1D(PrimR,fastest=fastestR)

LambdaMax = MAX(ABS(fastestL),ABS(fastestR))

Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolverByRusanov
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
