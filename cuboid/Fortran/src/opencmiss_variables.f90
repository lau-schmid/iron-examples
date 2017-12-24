
MODULE OPENCMISS_VARIABLES

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  LOGICAL :: CustomProfilingEnabled !< If custom profiling is compiled in
  LOGICAL :: TauProfilingEnabled !< If TAU profiling is compiled in

  INTEGER(CMISSIntg) :: Ftype,fibre_nr,NearestGP,InElement

  LOGICAL :: less_info,fast_twitch

  REAL(CMISSRP) :: TimeStart, TimeInitFinshed, TimeStretchSimFinished, TimeMainSimulationStart, TimeMainSimulationFinished
  REAL(CMISSSP), DIMENSION(2) :: DurationSystemUser     ! For receiving user and system time
  REAL(CMISSSP) :: DurationTotal
  REAL(CMISSRP) :: StimValuePerNode
  
  INTEGER(CMISSIntg) :: CustomSolverConvergenceReasonParabolic = 0
  INTEGER(CMISSIntg) :: CustomSolverConvergenceReasonNewton = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsParabolic = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsParabolicMin = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsParabolicMax = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsNewton = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsNewtonMin = 0
  INTEGER(CMISSIntg) :: CustomSolverNumberIterationsNewtonMax = 0
  INTEGER(CMISSIntg) :: MonodomainSolverId = 2
  INTEGER(CMISSIntg) :: MonodomainPreconditionerId = 1

  INTEGER(CMISSIntg) :: ODESolverId = 1

  INTEGER(CMISSIntg), DIMENSION(10000,100) :: MotorUnitFiringTimes
  INTEGER(CMISSIntg), ALLOCATABLE :: InnervationZoneOffset(:)
  INTEGER(CMISSIntg) :: JunctionNodeNo
  INTEGER(CMISSIntg) :: MotorUnitFires, MotorUnitRank

  INTEGER(CMISSIntg), ALLOCATABLE :: MuDistribution(:)
  REAL(CMISSRP) :: CustomTimingOdeSolver, CustomTimingParabolicSolver, CustomTimingFESolver, CustomTimingFileOutputUser, &
    & CustomTimingFileOutputSystem
  REAL(CMISSRP) :: CustomTimingOdeSolverPreLoad, CustomTimingParabolicSolverPreLoad, CustomTimingFESolverPreLoad, &
    & CustomTimingFileOutputUserPreLoad, CustomTimingFileOutputSystemPreLoad
  REAL(CMISSRP) :: TimingExportEMGUser = 0_CMISSRP
  REAL(CMISSRP) :: TimingExportEMGSystem = 0_CMISSRP

  !--------------------------------------------------------------------------------------------------------------------------------
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3

  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=NumberOfSpatialCoordinates
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussPoints=3

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensionsFE=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesFE=2

  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponentsFE=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MonodomainMeshComponentNumber=1

  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumberM=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberFE=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberFE=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberFE=4
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberFE=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberFE=6
  
  
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberM=2
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=NumberOfSpatialCoordinates

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumberM=3
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumberM=4
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariablesM=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponentsM=3 !Am, Cm, Conductiity   !(scalar, since 1D)

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumberM=6
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesM=3

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariablesFE=2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfComponentsFE=NumberOfSpatialCoordinates+1


  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentUserNumberM=9
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsM1=1
  INTEGER(CMISSIntg), PARAMETER :: FieldIndependentNumberOfComponentsM2=6 !5

  INTEGER(CMISSIntg), PARAMETER :: FieldEquationsSetUserNumberM=10

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberM=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetsUserNumberFE=2

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: SolverDAEIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverParabolicIndex=2
  INTEGER(CMISSIntg), PARAMETER :: SolverFEIndex=1

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopMonodomainNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopElasticityNumber=2

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: EquationsSetIndexM,EquationsSetIndexFE
  INTEGER(CMISSIntg) :: CellMLIndex
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,ComponentNumber,ElementDomain
  INTEGER(CMISSIntg) :: FibreNo, k
  INTEGER(CMISSIntg) :: NumberFiringFibres

  INTEGER(CMISSIntg), ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: shortenModelIndex!,shortenModelIndex2
  INTEGER(CMISSIntg) :: StimComponent

!  REAL(CMISSRP) :: YVALUE
  REAL(CMISSRP) :: VALUE

  INTEGER(CMISSIntg) :: Err, NodeIdx


  !CMISS variables

  TYPE(cmfe_BasisType) :: QuadraticBasis,LinearBasis,LinearBasisM
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(cmfe_CellMLType) :: CellMLEnvironment
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoopMain
  TYPE(cmfe_ControlLoopType) :: ControlLoopM,ControlLoopFE
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystemFE,CoordinateSystemM,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: DecompositionFE
  TYPE(cmfe_DecompositionType) :: DecompositionM
  TYPE(cmfe_EquationsType) :: EquationsFE
  TYPE(cmfe_EquationsType) :: EquationsM
  TYPE(cmfe_EquationsSetType) :: EquationsSetFE
  TYPE(cmfe_EquationsSetType) :: EquationsSetM
  TYPE(cmfe_FieldType) :: EquationsSetFieldFE
  TYPE(cmfe_FieldType) :: GeometricFieldFE
  TYPE(cmfe_FieldType) :: DependentFieldFE
  TYPE(cmfe_FieldType) :: IndependentFieldFE
  TYPE(cmfe_FieldType) :: MaterialFieldFE
  TYPE(cmfe_FieldType) :: EquationsSetFieldM
  TYPE(cmfe_FieldType) :: DependentFieldM
  TYPE(cmfe_FieldType) :: GeometricFieldM
  TYPE(cmfe_FieldType) :: IndependentFieldM
  TYPE(cmfe_FieldType) :: MaterialFieldM
  TYPE(cmfe_FieldType) :: FibreField
  TYPE(cmfe_FieldType) :: CellMLModelsField
  TYPE(cmfe_FieldType) :: CellMLStateField
  TYPE(cmfe_FieldType) :: CellMLIntermediateField
  TYPE(cmfe_FieldType) :: CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMeshFE
  TYPE(cmfe_MeshType) :: MeshFE
  TYPE(cmfe_MeshType) :: MeshM
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: RegionFE,RegionM,WorldRegion
  TYPE(cmfe_SolverType) :: SolverDAE,SolverParabolic
  TYPE(cmfe_SolverType) :: SolverFE,LinearSolverFE
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsM,SolverEquationsFE
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_MeshElementsType) :: QuadraticElements
  TYPE(cmfe_MeshElementsType) :: LinearElements
  TYPE(cmfe_MeshElementsType) :: ElementsM


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

END MODULE OPENCMISS_VARIABLES