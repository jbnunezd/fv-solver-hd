!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
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
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: KappaP1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

InitialCondition = 211

PrimRefState1 = (/1.4,0.0,0.0,1.0/)
PrimRefState2 = (/1.4,0.0,0.0,1.0/)
PrimRefState3 = (/1.4,0.0,0.0,1.0/)
PrimRefState4 = (/1.4,0.0,0.0,1.0/)

SELECT CASE(InitialCondition)
  CASE(200) ! Constant State
    TEnd    = 1.0
    Kappa   = 1.4
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE(211) ! Riemann Problem
    TEnd    = 1.0
    Kappa   = 1.4
    nElemsX = 600
    nElemsY = 600
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    PrimRefState1 = (/0.5313,0.0,0.0,0.4/)
    PrimRefState2 = (/1.0,0.7276,0.0,1.0/)
    PrimRefState3 = (/0.8,0.0,0.0,1.0/)
    PrimRefState4 = (/1.0,0.0,0.7276,1.0/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE(214) ! Double Mach Reflection Problem
    TEnd    = 0.25
    Kappa   = 1.4
    nElemsX = 2400
    nElemsY = 600
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/4.0,1.0/)
    PrimRefState1 = (/8.0,8.25*COS(PI/6.0),-8.25*SIN(PI/6.0),116.5/)
    PrimRefState2 = (/1.4,0.0,0.0,1.0/)
    PrimRefState3 = (/1.4,0.0,0.0,1.0/)
    PrimRefState4 = (/8.0,8.25*COS(PI/6.0),-8.25*SIN(PI/6.0),116.5/)
    BoundaryConditionsType = (/11,4,11,3/)
  CASE(215) ! Implosion Problem
    TEnd    = 2.5
    Kappa   = 1.4
    nElemsX = 600
    nElemsY = 600
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/0.3,0.3/)
    PrimRefState1 = (/1.0,0.0,0.0,1.0/)
    PrimRefState2 = (/0.125,0.0,0.0,0.14/)
    BoundaryConditionsType = (/5,5,5,5/)
  CASE(216) ! Shock-Vortex Interaction
    TEnd    = 0.2
    Kappa   = 1.4
    nElemsX = 1200
    nElemsY = 600
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/2.0,1.0/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

CFL      = 0.95
KappaM1  = Kappa-1.0
KappaP1  = Kappa+1.0
sKappaM1 = 1.0/KappaM1

Reconstruction    = 4
ReconstructionFix = 3

WhichOutput  = 1
nOutputFiles = 100

VarNameVisu(1) = 'Density'
VarNameVisu(2) = 'VelocityX'
VarNameVisu(3) = 'VelocityY'
VarNameVisu(4) = 'Pressure'

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
