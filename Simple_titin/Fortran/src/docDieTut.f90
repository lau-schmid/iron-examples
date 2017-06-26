! --------------------------------------------------------------------------------------------
! style idea: call "amount" numbers 'N_', user numbers 'UN_'
! ( and specific numbers 'k_' (i.e. 'for identification', similar to UN) )
!
! d: dependent. es: equationsset. f: fibre. g: geometric. i: independent. im: intermediate. 
! m: material. ms: models. p: parameters. s: state.
! Comp(s): component(s). CS: coordinate system.  Cml: CellML. Decomp: decomposition.
! Equs: equations. FE: Finte Elasticity. M: Monodomain. Par: parabolic. R: region.
! Sol: solver. Var(s): variable(s).
! --------------------------------------------------------------------------------------------
!                     *fortran - max number of signs per line: 132*
  
!> Main program
PROGRAM TITINEXAMPLE

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

! general program variables -----------------------------.
! i/o file status variable ------------------------------|
  INTEGER(CMISSIntg) :: stat                             !
! input file and path information -----------------------|
  CHARACTER(len=64), PARAMETER :: pathname = "input/"    !
!  CAUTION: File contains the complete micro physics:    |
  CHARACTER(len=64), PARAMETER :: filename = trim(pathname)// &    !
  & "Aliev_Panfilov_Razumova_Titin_2016_10_10_b.cellml"  !
! cpu time analysis variables ---------------------------|
! For receiving user and system time:                    |
  REAL :: elapsed(2)                                     !
! Total time (sum of prev. 2):                           |
  REAL :: total                                          !
! -------------------------------------------------------*

! improving readability a bit: -----------------------------------------------------------------------.
! ------------------------- field types --------------------------------------------------------------|
  INTEGER(CMISSIntg), PARAMETER ::  Geo_type      = CMFE_FIELD_GEOMETRIC_TYPE                  !   = 1|
  INTEGER(CMISSIntg), PARAMETER ::  F_type        = CMFE_FIELD_FIBRE_TYPE                      !   = 2|
  INTEGER(CMISSIntg), PARAMETER ::  Gen_type      = CMFE_FIELD_GENERAL_TYPE                    !   = 3|
  INTEGER(CMISSIntg), PARAMETER ::  Mat_type      = CMFE_FIELD_MATERIAL_TYPE                   !   = 4|
  INTEGER(CMISSIntg), PARAMETER ::  GeoGen_type   = CMFE_FIELD_GEOMETRIC_GENERAL_TYPE          !   = 5|
! ------------------------- variable types -----------------------------------------------------------|
  INTEGER(CMISSIntg), PARAMETER ::  U_type        = CMFE_FIELD_U_VARIABLE_TYPE                 !   = 1|         
  INTEGER(CMISSIntg), PARAMETER ::  V_type        = CMFE_FIELD_V_VARIABLE_TYPE                 !   = 5|
  INTEGER(CMISSIntg), PARAMETER ::  U1_type       = CMFE_FIELD_U1_VARIABLE_TYPE                !   = 9|
  INTEGER(CMISSIntg), PARAMETER ::  U2_type       = CMFE_FIELD_U2_VARIABLE_TYPE                !   =13|
  INTEGER(CMISSIntg), PARAMETER ::  DELUDELN_type = CMFE_FIELD_DELUDELN_VARIABLE_TYPE          !   = 2|
! ------------------------- kind of interpolation (more: GRID_POINT_BASED=4 DATA_POINT_BASED=6) ------|
  INTEGER(CMISSIntg), PARAMETER ::  C_IP          = CMFE_FIELD_CONSTANT_INTERPOLATION          !   = 1|
  INTEGER(CMISSIntg), PARAMETER ::  EB_IP         = CMFE_FIELD_ELEMENT_BASED_INTERPOLATION     !   = 2|
  INTEGER(CMISSIntg), PARAMETER ::  NB_IP         = CMFE_FIELD_NODE_BASED_INTERPOLATION        !   = 3|
  INTEGER(CMISSIntg), PARAMETER ::  GPB_IP        = CMFE_FIELD_GAUSS_POINT_BASED_INTERPOLATION !   = 5|
! ----------------------------------------------------------------------------------------------------|
  INTEGER(CMISSIntg), PARAMETER ::  SET_type      = CMFE_FIELD_VALUES_SET_TYPE                 !   = 1!
! ----------------------------------------------------------------------------------------------------*

! ===========================================================================================.
! ==========================================================================================.|
! ---------  M O D E L   D E S C R I P T I O N ,   O p e n C M I S S   S E T U P  ---------:||
! ==========================================================================================*|
! ===========================================================================================*

! space variables ---------------------------------------------------------------.
! all lengths in [cm]                                                            |
  REAL(CMISSRP), PARAMETER      :: LENGTH = 1.0_CMISSRP ! X-direction            |
  REAL(CMISSRP), PARAMETER      :: WIDTH  = 1.0_CMISSRP ! Y-direction            |
  REAL(CMISSRP), PARAMETER      :: HEIGHT = 1.0_CMISSRP ! Z-direction            |
!                                                                                |
! variables for the numerical space discretisation ----------------------------. |
! finite element discretisation (global muscle)                           -----| |
  INTEGER(CMISSIntg)            :: N_GlXElements, N_GlYElements, N_GlZElements ! |
  INTEGER(CMISSIntg)            :: N_ElementsFE                                ! |
!                                                                              | |
! fibre discretisation (monodomain)                                       -----| |
  INTEGER(CMISSIntg)            :: N_NodesInXi1, N_NodesInXi2                  ! |
  INTEGER(CMISSIntg)            :: N_NodesInXi3                                ! |
  INTEGER(CMISSIntg)            :: N_NodesPerFibre                             ! |
  INTEGER(CMISSIntg), PARAMETER :: N_FibresInSeries = 1                        ! |
  INTEGER(CMISSIntg)            :: N_NodesM                                    ! |
  INTEGER(CMISSIntg)            :: N_ElementsM                                 ! |
! -----------------------------------------------------------------------------* |
! -------------------------------------------------------------------------------*

! time variables ----------------------------------------------------------------.
! all times in [ms]. Time axis: [0 < TIME_STOP <= TIME_STOP_2]                   |
  REAL(CMISSRP), PARAMETER      :: time_default = 0.0001_CMISSRP                 !
  REAL(CMISSRP)                 :: time ! time variable                          |
  REAL(CMISSRP), PARAMETER      :: TIME_STOP = time_default! muscle fix until-/ stretch possible from- now.
  REAL(CMISSRP), PARAMETER      :: TIME_STOP_2 = time_default! end of simulation |
! time span how long the muscle is stimulated. typically about 1/10 ms:          |
  REAL(CMISSRP), PARAMETER      :: STIM_STOP = time_default                      !
                                                                                 !
! numerical time step sizes ------------------------------------------.          |
  REAL(CMISSRP), PARAMETER      :: TIME_STEP           = time_default !          |
  REAL(CMISSRP), PARAMETER      :: ODE_TIME_STEP       = time_default !          |
  REAL(CMISSRP), PARAMETER      :: PDE_TIME_STEP       = time_default !          |
  REAL(CMISSRP), PARAMETER      :: ELASTICITY_TIME_STEP= time_default !          |
! --------------------------------------------------------------------*          |
! -------------------------------------------------------------------------------*

! numerical constant -------------------------------------------------.
  REAL(CMISSRP), PARAMETER      :: tol = 1.0E-8_CMISSRP               !
! --------------------------------------------------------------------*

! --------------------------------------------------------------------------------------------
  INTEGER(CMISSIntg), PARAMETER :: OUTPUT_FREQUENCY = 1 !1
! --------------------------------------------------------------------------------------------  
  REAL(CMISSRP)                 :: stretch_sarcolength_ratio = 1.0_CMISSRP
! physical parameters ------------------------------------------------------.
! to set the mechanical problem trivial, i.e stresses are ruled out.  ------|
  REAL(CMISSRP), PARAMETER      :: P_max = 0.0_CMISSRP ! [N/cm^2]           |
! 1: With Actin-Titin Interaction 0: No Actin-Titin Interactions      ------|
  REAL(CMISSRP), PARAMETER      :: TK_lin_param = 0.0_CMISSRP               !
! stimulation current                                                 ------|
  REAL(CMISSRP)                 :: STIM_VALUE                               !
! condctivity in [mS/cm]                                              ------|
  REAL(CMISSRP), PARAMETER      :: CONDUCTIVITY = 0.5_CMISSRP               !
! surface area to volume ratio in [cm^-1]                             ------|
  REAL(CMISSRP), PARAMETER      :: Am           = 1.0_CMISSRP               !
! membrane capacitance in [uF/cm^2]                                   ------|
  REAL(CMISSRP), PARAMETER      :: Cm_fast      = 1.0_CMISSRP               !
! maximum contraction velocity in [cm/ms]                             ------|
  REAL(CMISSRP), PARAMETER      :: Vmax         = -0.015_CMISSRP ! [m/s]    |
! material constants, CAUTION - what are the units???                 ------|
  REAL(CMISSRP), PARAMETER      :: alpha = 5.0_CMISSRP                      !
  REAL(CMISSRP), PARAMETER      :: beta  = 0.01_CMISSRP                     !
  REAL(CMISSRP), PARAMETER      :: gama  = 1.0_CMISSRP                      !
! c_{10}, c_{01} parameters [N/cm^2] as from Heidlauf, Table 6.1::
  REAL(CMISSRP), PARAMETER :: MooneyRivlin1st = 0.0000000000635201_CMISSRP  !
  REAL(CMISSRP), PARAMETER :: MooneyRivlin2nd = 0.3626712895523322_CMISSRP  !
  REAL(CMISSRP), PARAMETER, DIMENSION(4) :: MAT_FE = & ![N/cm^2]      ------|
  & [alpha*MooneyRivlin1st, alpha*MooneyRivlin2nd, &                        !
  &  beta*1.074519943356914_CMISSRP, gama*9.173311371574769_CMISSRP]        !
! initial condition:?                                                 ------|
  REAL(CMISSRP) :: INIT_PRESSURE = -2.0_CMISSRP*MAT_FE(2)-MAT_FE(1)         !
! --------------------------------------------------------------------------*
 
  
! ----- internal numberings found in 8 structs:-----------------------------------------------
! CoordinateSystem (Coordinate Systems),
! Region (Regions),
! *Basis (Basis Functions),
! GeneratedMeshUser,
! Mesh,
! Decomposition,
! *Field*,
! EquationsSets(!!s, so not: EquationsSet), Problem (Problems),
! --------------------------------------------------------------------------------------------
! "missing" from 'top_structure.png': (Base System), (Computational Environment)
! --------------------------------------------------------------------------------------------
! NOT found in "Solver" or "ControlLoop"
! --------------------------------------------------------------------------------------------
! Contained by 'CSWorld' (w or w/o UserNumber!?)::
  INTEGER(CMISSIntg), PARAMETER :: UN_CSFE = 1 ! was 'CoordinateSystemUserNumberFE'
  INTEGER(CMISSIntg), PARAMETER :: UN_CSM  = 2 ! was 'CoordinateSystemUserNumberM'
  INTEGER(CMISSIntg), PARAMETER :: N_Dim   = 3 ! to work in 3D. was 'NumberOfSpatialCoordinates'
  INTEGER(CMISSIntg), PARAMETER :: k_CSFEDim = 3
  INTEGER(CMISSIntg), PARAMETER :: k_CSMDim = 1
! Contained by 'WorldRegion' (w or w/o UserNumber!?)::
  INTEGER(CMISSIntg), PARAMETER :: UN_RFE  = 1  ! was 'RegionUserNumberFE'
  INTEGER(CMISSIntg), PARAMETER :: UN_RM   = 2

  INTEGER(CMISSIntg), PARAMETER :: UN_QuadBasisFE = 1
  INTEGER(CMISSIntg), PARAMETER :: UN_LinBasisFE  = 2
  INTEGER(CMISSIntg), PARAMETER :: UN_LinBasisM   = 3
  
  ! replaced (2 items) by N_Dim
  INTEGER(CMISSIntg), PARAMETER :: N_GaussPointsFE = 3 ! in context of a quadratic basis (and a 'common' linear one) 
  INTEGER(CMISSIntg), PARAMETER :: N_GaussPointsM = 2 ! in context of a linear basis 

  INTEGER(CMISSIntg), PARAMETER :: UN_GenMesh = 1        ! 1--.  
  INTEGER(CMISSIntg), PARAMETER :: UN_MeshFE = 1         !    '-1
  INTEGER(CMISSIntg), PARAMETER :: UN_MeshM  = 2         !      2
  INTEGER(CMISSIntg), PARAMETER :: N_MeshMDim = k_CSMDim

! not needed, since set automatically via the Generated Mesh:: 
!  INTEGER(CMISSIntg), PARAMETER :: N_MeshFEComps = 2 ! war 'NumberOfMeshComponentsFE'
!  INTEGER(CMISSIntg), PARAMETER :: N_MeshFEDim  = k_CSFEDim
  INTEGER(CMISSIntg), PARAMETER :: N_MeshMComps = 1 ! neue Variable. ersetzt eine '1'
! ! usage in subroutine 'cmfe_Field_ComponentMeshComponentSet(.,.,., *Comp* , Err)' ! ! :::
  INTEGER(CMISSIntg), PARAMETER :: MeshFEComp1 = 1 ! war 'QuadraticMeshComponentNumber'
  INTEGER(CMISSIntg), PARAMETER :: MeshFEComp2 = 2 ! war 'LinearMeshComponentNumber'
  INTEGER(CMISSIntg), PARAMETER :: MeshMComp1 = 1  ! war 'MonodomainMeshComponentNumber'
  
  INTEGER(CMISSIntg), PARAMETER :: UN_DecompFE = 1
  INTEGER(CMISSIntg), PARAMETER :: UN_DecompM  = 2

  INTEGER(CMISSIntg), PARAMETER :: UN_gFieldFE          = 1 !### (### are 'Field user numbers')
  INTEGER(CMISSIntg), PARAMETER :: Type_gFieldFE = Geo_type
  INTEGER(CMISSIntg), PARAMETER :: N_gFieldFEVars  = 1
  INTEGER(CMISSIntg), PARAMETER :: N_gFieldFEComps = N_Dim
  INTEGER(CMISSIntg), PARAMETER :: Type_gFieldFEVar = U_type
  
  INTEGER(CMISSIntg), PARAMETER :: UN_gFieldM           = 2 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_gFieldM = Geo_type
  INTEGER(CMISSIntg), PARAMETER :: N_gFieldMVars = 1  
  INTEGER(CMISSIntg), PARAMETER :: N_gFieldMComps = 1 
  INTEGER(CMISSIntg), PARAMETER :: Type_gFieldMVar = U_type

  INTEGER(CMISSIntg), PARAMETER :: UN_fField            = 3 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_fField = F_type
  INTEGER(CMISSIntg), PARAMETER :: N_fFieldVars = 1
  INTEGER(CMISSIntg), PARAMETER :: N_fFieldComps = 3 
  INTEGER(CMISSIntg), PARAMETER :: Type_fFieldVar = U_type

  INTEGER(CMISSIntg), PARAMETER :: UN_mFieldM           = 4 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_mFieldM = Mat_type
  INTEGER(CMISSIntg), PARAMETER :: N_mFieldMVars = 1
  INTEGER(CMISSIntg), PARAMETER :: N_mFieldMComps = 3 !Am, Cm, Conductivity   !(scalar, since 1D)
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_mFieldMComps)      :: IP_mFieldMVar = [C_IP, NB_IP, C_IP]
  INTEGER(CMISSIntg), PARAMETER :: Type_mFieldMVar = U_type 

  INTEGER(CMISSIntg), PARAMETER :: UN_mFieldFE          = 5 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_mFieldFE = Mat_type
  INTEGER(CMISSIntg), PARAMETER :: N_mFieldFEVars = 2
  INTEGER(CMISSIntg), PARAMETER :: N_mFieldFEVar1Comps = 6 ! (material properties) was '5'
  INTEGER(CMISSIntg), PARAMETER :: N_mFieldFEVar2Comps = 1 ! (gravity?)
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_mFieldFEVar1Comps) :: IP_mFieldFEVar1 = [C_IP, C_IP, C_IP, C_IP, C_IP, C_IP] ! not automated!
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_mFieldFEVar2Comps) :: IP_mFieldFEVar2 = [C_IP] ! not automated!
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_mFieldFEVars) :: Type_mFieldFEVar = [U_type, V_type]
  
  INTEGER(CMISSIntg), PARAMETER :: UN_dFieldM           = 6 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_dFieldM = Gen_type
  INTEGER(CMISSIntg), PARAMETER :: N_dFieldMVars=3
  INTEGER(CMISSIntg), PARAMETER :: N_dFieldMVar1Comps=1
  INTEGER(CMISSIntg), PARAMETER :: N_dFieldMVar2Comps=1
  INTEGER(CMISSIntg), PARAMETER :: N_dFieldMVar3Comps=3
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_dFieldMVars) :: Type_dFieldMVar = [U_type, DELUDELN_type, V_type]

  INTEGER(CMISSIntg), PARAMETER :: UN_dFieldFE          = 7 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_dFieldFE = GeoGen_type
  INTEGER(CMISSIntg), PARAMETER :: N_dFieldFEVars=2            ! (state and derivative of the state?)
  INTEGER(CMISSIntg), PARAMETER :: N_dFieldFEComps = N_Dim + 1 ! (displacement (3D -> 3 components), pressure (1 component)) each variable with four components? 
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_dFieldFEVars) :: Type_dFieldFEVar = [U_type, DELUDELN_type] ! why not DELUDELT?

  INTEGER(CMISSIntg), PARAMETER :: UN_iFieldFE          = 8 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_iFieldFE = Gen_type
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldFEVars=2
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldFEVar1Comps=6 ! 1:? 2:titin force (unbound) 3:titin force (bound) 4:titin force in XF-direction (unbound) 5:titin force in XF-direction (bound) 6:activation for titin
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldFEVar2Comps=4
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldFEVar1Comps) :: IP_iFieldFEVar1 = [GPB_IP,GPB_IP,GPB_IP,GPB_IP,GPB_IP,GPB_IP] ! not automated!
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldFEVar2Comps) :: IP_iFieldFEVar2 = [EB_IP, EB_IP, EB_IP, EB_IP] ! not automated!
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldFEVars) :: Type_iFieldFEVar = [U_type, V_type]
  
  INTEGER(CMISSIntg), PARAMETER :: UN_iFieldM           = 9 !###
  INTEGER(CMISSIntg), PARAMETER :: Type_iFieldM = Gen_type
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldMVars=4
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldMVar1Comps=6 !1:normalised sarcomere-based active stress 2:for the unbound titin stress 3:for the bound titin stress 4:for the unbound titin stress in the XF-direction 5:for the bound titin stress in the XF-direction 6:for the titin activation (stress without force-length relation)
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldMVar2Comps=5 !1:motor unit number 2:fibre type 3:fibre number 4:nearest Gauss point 5:in element number (LOCAL NODE NUMBERING!!!)
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldMVar3Comps=4 !1:sarcomere half length 2:inital sarcomere half length 3:initial node distance 4:sarcomere half length when activation starts
  INTEGER(CMISSIntg), PARAMETER :: N_iFieldMVar4Comps=6 !1:old node distance 2:maximum contraction velocity 3:relative contraction velocity 4:TK old node distance 2 5:TK old node distance 3 6:TK old node distance 4
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldMVar1Comps) :: IP_iFieldMVar1 = [NB_IP, NB_IP, NB_IP, NB_IP, NB_IP, NB_IP]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldMVar2Comps) :: IP_iFieldMVar2 = [NB_IP, NB_IP, NB_IP, NB_IP, NB_IP]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldMVar3Comps) :: IP_iFieldMVar3 = [NB_IP , C_IP, C_IP, NB_IP]
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldMVar4Comps) :: IP_iFieldMVar4 = [NB_IP, C_IP, NB_IP, NB_IP, NB_IP, NB_IP] 
  INTEGER(CMISSIntg), PARAMETER, DIMENSION(N_iFieldMVars) :: Type_iFieldMVar = [U_type, V_type, U1_type, U2_type]
  
  INTEGER(CMISSIntg), PARAMETER :: UN_esFieldM      = 10 !###
  INTEGER(CMISSIntg), PARAMETER :: UN_esFieldFE     = 11 !###

  INTEGER(CMISSIntg), PARAMETER :: UN_Cml = 1 ! UN of struct, which is needed, to set up the CellML fields
  INTEGER(CMISSIntg), PARAMETER :: UN_msFieldCml        = 12 !###
  INTEGER(CMISSIntg), PARAMETER :: UN_sFieldCml         = 13 !###
  INTEGER(CMISSIntg), PARAMETER :: UN_imFieldCml        = 14 !###
  INTEGER(CMISSIntg), PARAMETER :: UN_pFieldCml         = 15 !###

  INTEGER(CMISSIntg), PARAMETER :: UN_EquSetM  = 1
  INTEGER(CMISSIntg), PARAMETER :: UN_EquSetFE = 2

  INTEGER(CMISSIntg), PARAMETER :: UN_Problem=1

  INTEGER(CMISSIntg), PARAMETER :: idx_SolDAE = 1 ! does '1' mean DAE or DAE-version or any of these?!
  INTEGER(CMISSIntg), PARAMETER :: idx_SolPar = 2
  INTEGER(CMISSIntg), PARAMETER :: idx_SolFE  = 1

  INTEGER(CMISSIntg), PARAMETER :: k_ControlLoopM=1
  INTEGER(CMISSIntg), PARAMETER :: k_ControlLoopFE=2
  
! Program types
  
! Program variables

  INTEGER(CMISSIntg) :: indx_EqusSetM,idx_EqusSetFE ! set by program as needed by program
  INTEGER(CMISSIntg) :: idx_Cml ! set by program as needed by program
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: N_CNodes,N_Doms,k_CNode
  
  INTEGER(CMISSIntg) :: k_Node,NodeDomain,node_idx,domain_idx,ElementDomain !,ComponentNumber

  INTEGER(CMISSIntg) :: i,j,k,elem_idx,node,iter,comp_idx !,NumberOfElementsPerElasticityElement,my_node_idx,node2

  INTEGER(CMISSIntg), ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: idx_CmlSH!,shortenModelIndex2
  INTEGER(CMISSIntg) :: stimcomponent

  REAL(CMISSRP) :: VALUE
  
  INTEGER(CMISSIntg) :: Err

  LOGICAL                       :: independent_field_auto_create = .FALSE.

  !CMISS variables ---------------------------------------------------------------------------
  !  note: (most) objects of cmfe_*type contain (exactly) a private pointer to an object
  !  of their original type. I.e.: (A is of cmfe_BooType) ==> (A%boo is of Boo_Type).
  !CMISS variables
  TYPE(cmfe_BasisType)          :: QuadBasisFE,LinBasisFE,LinBasisM
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsM,BoundaryConditionsFE
  TYPE(cmfe_CellMLType)         :: Cml
  TYPE(cmfe_CellMLEquationsType) :: EqusCml
  TYPE(cmfe_ControlLoopType)    :: ControlLoopMain
  TYPE(cmfe_ControlLoopType)    :: ControlLoopM,ControlLoopFE
  TYPE(cmfe_CoordinateSystemType) :: CSFE,CSM,CSWorld
  TYPE(cmfe_DecompositionType)  :: DecompFE,DecompM
  TYPE(cmfe_EquationsType)      :: EqusM,EqusFE
  TYPE(cmfe_EquationsSetType)   :: EqusSetM,EqusSetFE
  TYPE(cmfe_FieldType)          :: esFieldM, gFieldM, dFieldM, iFieldM, mFieldM
  TYPE(cmfe_FieldType)          :: esFieldFE,gFieldFE,dFieldFE,iFieldFE,mFieldFE
  TYPE(cmfe_FieldType)          :: fField
  TYPE(cmfe_FieldType)          :: msFieldCml,sFieldCml,imFieldCml,pFieldCml
  TYPE(cmfe_FieldsType)         :: Fields ! (for some kind of data export)
  TYPE(cmfe_GeneratedMeshType)  :: GenMesh
  TYPE(cmfe_MeshType)           :: MeshFE,MeshM
  TYPE(cmfe_ProblemType)        :: Problem
  TYPE(cmfe_RegionType)         :: RFE,RM,RWorld !WldRg. is parent_region of other two.
  TYPE(cmfe_SolverType)         :: SolverDAE,SolverParabolic
  TYPE(cmfe_SolverType)         :: SolverFE,LinearSolverFE
  TYPE(cmfe_SolverEquationsType) :: SolEqusM,SolEqusFE
  TYPE(cmfe_NodesType)          :: NodesM
  ! 'QuadraticElms, LinearElms' not needed since mesh created automatically via GeneratedMesh
  !TYPE(cmfe_MeshElementsType) :: QuadraticElements
  !TYPE(cmfe_MeshElementsType) :: LinearElements
  TYPE(cmfe_MeshElementsType)   :: ElementsM

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
!---------------------------------------------------------------------------------------------
!---  E N D   O F   V A R I A B L E   D E C L A R A T I O N  ---------------------------------
!---  N O W   M O R E   C O M P L E X   S T R U C T S   A R E   C R E A T E D  --------||-----
!---------------------------------------------------------------------------------------*

!#############################################################################################
!################################### SET DISCRETISATION STUFF ################################
!#############################################################################################
! variables for the numerical space discretisation ----------------------------. |############
! finite element discretisation (global muscle)                           -----| |############
  N_GlXElements = 1                                                            ! |############
  N_GlYElements = 1                                                            ! |############
  N_GlZElements = 1                                                            ! |############
  N_ElementsFE = N_GlXElements*N_GlYElements*N_GlZElements                     ! |############
!#############################################################################################
! fibre discretisation (monodomain)                                       -----| |############
  N_NodesInXi1 = 2                                                       ! |############
  N_NodesInXi2 = 1                                                       ! |############
  N_NodesInXi3 = 1                                                       ! |############
  N_NodesPerFibre = (N_NodesInXi1-1)*N_GlXElements/N_FibresInSeries+1    ! |############
  N_NodesM = N_NodesPerFibre*N_FibresInSeries*N_GlYElements*N_GlZElements* &   ! |############
           & N_NodesInXi2*N_NodesInXi3                             ! |############
  N_ElementsM = (N_NodesPerFibre-1)*N_FibresInSeries*N_GlYElements* &          ! |############
              & N_GlZElements*N_NodesInXi2*N_NodesInXi3            ! |############
! -----------------------------------------------------------------------------* |############
! -------------------------------------------------------------------------------*############
!#############################################################################################
! input file information --------------------------------.                        ############
! pathname = "input/"                                     !                        ############
! CAUTION: File contains the complete micro physics:     |                        ############
! filename = trim(pathname)//"Aliev_Panfilov_Razumova_Titin_2016_10_10_b.cellml"  !############
!--------------------------------------------------------*                        ############
!#############################################################################################
!################################### SET DISCRETISATION STUFF ################################
!#############################################################################################


!======================================================================================||
!============================== Build up the cmfe structs =============================||
!======================================================================================||
! ----------------------------------------------------------------------------------|
! -Intialise OpenCMISS--------------------------------------------------------------|
  CALL cmfe_Initialise(CSWorld,RWorld,Err)
! -Trap errors----------------------------------------------------------------------|
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
! - MPI information. Get computational nodes information: --------------------------|
! - first: 'Number of CNodes'                                                       |
  CALL cmfe_ComputationalNumberOfNodesGet(N_CNodes,Err)
! - second: 'the specific number assigned to THIS CNode'                            |
  CALL cmfe_ComputationalNodeNumberGet(k_CNode,Err)
! - some (I/)O setting -------------------------------------------------------------|
  CALL cmfe_OutputSetOn("s_titin",Err)
! ----------------------------------------------------------------------------------|
! the domain decomposition at this moment: for each CNode a separate domain:        |
  N_Doms = N_CNodes ! CARE! >1 might throw errors!                                  |
! ----------------------------------------------------------------------------------|
! = MPI BROADCAST ==================================================================|
! Broadcast the number of elements in the X,Y,Z directions and the number of        |
! domain decompositions to the other computational nodes                            |
  CALL MPI_BCAST(N_GlXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(N_GlYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(N_GlZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(N_Doms,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
! ----------------------------------------------------------------------------------|


! COORDINATE SYSTEM:===================================================================||
! ----------------------------------------------------------------------------------|
! --------BEGIN---------------------------------------------------------------------|
! ------- Level2 - coordinate systems ------ note: no correlation between them! ----|
! -------   contained by WorldCoordinateSystem CSWorld -----------------------------|
! ---------- 1: --------------------------------------------------------------------|
! Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CSFE,Err)
  CALL cmfe_CoordinateSystem_CreateStart(UN_CSFE,CSFE,Err)
! Set the coordinate system to be 3D - was '3' instead of k_CSFEDim
  CALL cmfe_CoordinateSystem_DimensionSet(CSFE,k_CSFEDim,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CSFE,Err)
! ---------- 2: --------------------------------------------------------------------|
! Create a 1D coordinate System
  CALL cmfe_CoordinateSystem_Initialise(CSM,Err)
  CALL cmfe_CoordinateSystem_CreateStart(UN_CSM,CSM,Err)
! Set the coordinate system to be 1D - was '1' instead of k_CSMDim
  CALL cmfe_CoordinateSystem_DimensionSet(CSM,k_CSMDim,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CSM,Err)
! --------END-----------------------------------------------------------------------|
! ----------------------------------------------------------------------------------|
  
! REGION:==============================================================================||
! ----------------------------------------------------------------------------------|
! --------BEGIN---------------------------------------------------------------------|
! ------- Level2 - regions ---------------------------------------------------------|
! --------  contained by WorldRegion RWorld ----------------------------------------|
! ----------1: ---------------------------------------------------------------------|
  CALL cmfe_Region_Initialise(RFE,Err)
  CALL cmfe_Region_CreateStart(UN_RFE,RWorld,RFE,Err)
  CALL cmfe_Region_CoordinateSystemSet(RFE,CSFE,Err)
  CALL cmfe_Region_LabelSet(RFE,"Region3D",Err)
  CALL cmfe_Region_CreateFinish(RFE,Err)
! -------- 2: ----------------------------------------------------------------------|
  CALL cmfe_Region_Initialise(RM,Err)
  CALL cmfe_Region_CreateStart(UN_RM,RWorld,RM,Err)
  CALL cmfe_Region_CoordinateSystemSet(RM,CSM,Err)
  CALL cmfe_Region_LabelSet(RM,"Region1D",Err)
  CALL cmfe_Region_CreateFinish(RM,Err)
! --------END-----------------------------------------------------------------------|
! ----------------------------------------------------------------------------------|

! BASIS FUNCTIONS:=====================================================================||
  ! --------------------------------------------------------------------------------|
  ! -----BEGIN ---------------------------------------------------------------------|
  ! ----- Level1 - basis functions--------------------------------------------------|
  != CAUTION =======================================================================|
  ! QuadBasisFE and LinBasisFE need to: 1. be of same type;  2. have the same amount|
  ! of xi. To be able to use SUBROUTINE cmfe_GeneratedMesh_BasisSet(), later.       |
  ! ------- 1: ---------------------------------------------------------------------|
  ! Define basis 1 - tri-Quadratic Lagrange, FOR THE 3D GRID
  CALL cmfe_Basis_Initialise(QuadBasisFE,Err)
  CALL cmfe_Basis_CreateStart(UN_QuadBasisFE,QuadBasisFE,Err)
  CALL cmfe_Basis_TypeSet(QuadBasisFE,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadBasisFE,k_CSFEDim,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadBasisFE,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadBasisFE, &
   & [N_GaussPointsFE,N_GaussPointsFE,N_GaussPointsFE],Err)
  CALL cmfe_Basis_CreateFinish(QuadBasisFE,Err)
  ! ------ 2: ----------------------------------------------------------------------|
  ! Define basis 2 - tri-Linear Lagrange, FOR THE 3D GRID
  CALL cmfe_Basis_Initialise(LinBasisFE,Err)
  CALL cmfe_Basis_CreateStart(UN_LinBasisFE,LinBasisFE,Err)
  CALL cmfe_Basis_TypeSet(LinBasisFE,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinBasisFE,k_CSFEDim,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinBasisFE,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
   & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinBasisFE, &
   & [N_GaussPointsFE,N_GaussPointsFE,N_GaussPointsFE],Err)
  CALL cmfe_Basis_CreateFinish(LinBasisFE,Err)
  ! ------ 3: ----------------------------------------------------------------------|
  ! Define basis 3 - tri-Linear Lagrange, FOR THE 1D GRID
  CALL cmfe_Basis_Initialise(LinBasisM,Err)
  CALL cmfe_Basis_CreateStart(UN_LinBasisM,LinBasisM,Err)
  CALL cmfe_Basis_TypeSet(LinBasisM,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinBasisM,k_CSMDim,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinBasisM,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinBasisM,[N_GaussPointsM],Err)
  CALL cmfe_Basis_CreateFinish(LinBasisM,Err)
  ! ------END-----------------------------------------------------------------------|
  ! --------------------------------------------------------------------------------|

! MESHES:==============================================================================||
  ! --------------------------------------------------------------------------------|
  ! -----BEGIN ---------------------------------------------------------------------|
  ! NOTE: 'Mesh'es as well as 'GeneratedMesh'es are part of a region struct --------|
  ! --------------------------------------------------------------------------------|
  ! Create a Level1-'mesh' via a "containing" Level1-'generated mesh'.
  !  / I think both meshes are level 1 and there is an other struct 'gen.mesh' which is \
  ! (  level 1 as well, that is used to create the mesh struct more easily. while it is  )
  !  \ possibl to have a mesh embedding others we don't use this here...                /
  !---------------------------------------------------------------------------------|
  ! Mesh will contain N elements, N = N_GlXElements x NGYE x NGZE.
  
! GENERATEDMESHES:====================================||
  ! -----INTERMEDIATE BEGIN-------------------------|
  ! ---- Level1-GeneratedMesh ----------------------|
  ! ------ 1: --------------------------------------|
  ! Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GenMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(UN_GenMesh,RFE,GenMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GenMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GenMesh,[QuadBasisFE,LinBasisFE],Err)
  !-------------------------------------------------|
  ! Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GenMesh, [LENGTH,WIDTH,HEIGHT], Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GenMesh, &
    & [N_GlXElements, N_GlYElements, N_GlZElements], Err)
  ! -----INTERMEDIATE INTERRUPT---------------------|
  
! ---- mesh 1: ---------------------------------------------------------------------|
! create a 'Mesh'
  CALL cmfe_Mesh_Initialise(MeshFE, Err)
! ----------------------------------------------------------------------------------|
  
  ! -----INTERMEDIATE CONTINUE----------------------|
  ! Finish the creation of 'the generated mesh' in
  ! the region and initialise all MeshFE mesh parts
  ! automatically:
  CALL cmfe_GeneratedMesh_CreateFinish(GenMesh,UN_MeshFE,MeshFE,Err)
  ! here, NumberOfComponents/ -Elements is set (N_comps="SIZE(REGULAR_MESH%BASES)") 
  ! as well as the number of dimensions as "GenMesh%Basis%NUMBER_OF_XI
  ! -----INTERMEDIATE END---------------------------|
!=====================================================||

! ---- mesh 2: ---------------------------------------------------------------------|  
! Create a mesh in the region directly
  CALL cmfe_Mesh_Initialise(MeshM,Err)
  CALL cmfe_Mesh_CreateStart(UN_MeshM,RM,N_MeshMDim,MeshM,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(MeshM,N_MeshMComps,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(MeshM,N_ElementsM,Err)

  CALL cmfe_MeshElements_Initialise(ElementsM,Err)
  CALL cmfe_MeshElements_CreateStart(MeshM,MeshMComp1,LinBasisM,ElementsM,Err)

! Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(NodesM,Err)
  CALL cmfe_Nodes_CreateStart(RM,N_NodesM,NodesM,Err)
  CALL cmfe_Nodes_CreateFinish(NodesM,Err)
! ----------------------------------------------------------------------------------|
  elem_idx=0
  DO node=1,N_NodesM
    IF(mod(node,N_NodesPerFibre)==0) CYCLE ! treat last node in interval different
    elem_idx=elem_idx+1
    CALL cmfe_MeshElements_NodesSet(ElementsM,elem_idx,[node,node+1],Err) ! todo potential speedup
  ENDDO    
! ----------------------------------------------------------------------------------|
! this implements 1D meshes: node_i--(elem_idx)--node_{i+1} ------------------------|
! ----------------------------------------------------------------------------------|
! example mesh for: N_NodesPerFibre = 5, N_NodesM = 15 -----------------|
!
!     [1---(1)--2---(2)---3---(3)---4---(4)---5 ], mesh part 1 of MeshM        -----|
!
!     [6---(5)--7---(6)---8---(7)---9---(8)---10], mesh part 2 of MeshM        -----|
!
!     [11--(9)--12--(10)--13--(11)--14--(12)--15], mesh part 3 of MeshM        -----|
! ----------------------------------------------------------------------------------|
  write(*,*) "Finished setting up 1D elements"
  CALL cmfe_MeshElements_CreateFinish(ElementsM,Err)
  CALL cmfe_Mesh_CreateFinish(MeshM,Err)
! --------END-----------------------------------------------------------------------|
! ----------------------------------------------------------------------------------|
  
! DECOMPOSITIONS:======================================================================||
! ----------------------------------------------------------------------------------|
! -------BEGIN ---------------------------------------------------------------------|
! --- Note: Decompositions are assigned to exactly one mesh. A decomposition tells -|
! --- whether and how a mesh is devided into (sub)domains. It maps every element of |
! --- a mesh to exactly one such domain. A mesh could have multiple decompositions. |
! ----Level1 - decompostions--------------------------------------------------------|
! ------ 1: ------------------------------------------------------------------------|
! Create a decomposition (for the MeshFE mesh)
  CALL cmfe_Decomposition_Initialise(DecompFE,Err)
! Link it to the the MeshFE mesh ---------------------------------------------------|
  CALL cmfe_Decomposition_CreateStart(UN_DecompFE,MeshFE,DecompFE,Err)
! spread the elements "in a fair manner" over the domains --------------------------|
  IF(N_Doms>1) THEN           ! Actual decomposition CASE *****
    CALL cmfe_Decomposition_TypeSet(DecompFE,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    elem_idx=0                      ! e=0         (e.g.: N_Doms=4, N_ElementsFE=16)
    DO domain_idx=0,N_Doms-1        ! d=0      |1      |2         |3
      DO iter=1,N_ElementsFE/N_Doms ! i=1,2,3,4|1,2,3,4|1, 2, 3, 4| 1, 2, 3, 4
        elem_idx=elem_idx+1         ! e=1,2,3,4|5,6,7,8|9,10,11,12|13,14,15,16
        ! assign a domain to each element: (e.g. elements 5 to 8 are in domain 1)
        CALL cmfe_Decomposition_ElementDomainSet(DecompFE,elem_idx,domain_idx,Err)
      ENDDO
    ENDDO
  ELSE                        ! Trivial decomposition CASE *****
    CALL cmfe_Decomposition_TypeSet(DecompFE,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err) 
  ENDIF
  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompFE,N_Doms,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(DecompFE,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(DecompFE,Err)
! ----------------------------------------------------------------------------------|
! ------ 2: ------------------------------------------------------------------------|
! Create a second decomposition (one for the MeshM mesh)
  CALL cmfe_Decomposition_Initialise(DecompM,Err)
! Link it to the MeshM mesh --------------------------------------------------------|
  CALL cmfe_Decomposition_CreateStart(UN_DecompM,MeshM,DecompM,Err)

  IF(N_Doms>1) THEN          ! Actual decomposition CASE *****
! #################################################################################
! ######################### S P E C I A L   C A S E ! #############################
! ######## this only works for a domain decomposition of "2x2x2" domains  #########
! #################################################################################
  ! Create a domain for every computational node, map the mesh elements respectively
    CALL cmfe_Decomposition_TypeSet(DecompM,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    elem_idx=0
    DO domain_idx=0,N_Doms-1            !  0         1         2         3
      DO iter=1,N_ElementsFE/N_Doms/2   !  1    2    1    2    1    2    1    2    what's the sense in this index?? das funktioniert so eventuell, wenn in alle richtungen 2 domains sind
        DO i=1,N_NodesInXi3       !  1 2  1 2  1 2  1 2  1 2  1 2  1 2  1 2
          DO j=1,N_NodesInXi2     !  1212 1212 1212 1212 1212 1212 1212 1212        
            DO k=1,N_NodesPerFibre-1    ! 
              elem_idx=elem_idx+1       !  1234 5678 9abc defg hijk lmno pqrs tuvw 
              CALL cmfe_Decomposition_ElementDomainSet(DecompM,elem_idx,domain_idx,Err)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! Check whether all elements have been assigned to a domain:
    IF(elem_idx/=N_ElementsM) THEN
      WRITE(*,*) "Error in setting up the decomposition for monodomain!"
      STOP
    ENDIF
    ! ------domains created---------------------------------------------------------|
  ELSE                       ! Trivial decomposition CASE *****
    CALL cmfe_Decomposition_TypeSet(DecompM,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  ENDIF
  CALL cmfe_Decomposition_NumberOfDomainsSet(DecompM,N_Doms,Err)
  CALL cmfe_Decomposition_CreateFinish(DecompM,Err)
! -----END--------------------------------------------------------------------------|
! ----------------------------------------------------------------------------------|

! FIELDS: =============================================================================||
! ----------------------------------------------------------------------------------|
! -----BEGIN -----------------------------------------------------------------------|
! !=================================================================================|
! !  F I N I T E   E L A S T C I T Y                                                |
! !=================================================================================|
! there are 4 fields necessary: geometric, material, dependent and independent.-----|
! ----------------------------------------------------------------------------------|
! Create a geometric field for finite elasticity - quadratic interpolation
! (is this "create function, set domain/codomain and assign decomposition"?) -------?
  CALL cmfe_Field_Initialise(gFieldFE,Err)
  CALL cmfe_Field_CreateStart(UN_gFieldFE,RFE,gFieldFE,Err)
  CALL cmfe_Field_TypeSet(gFieldFE,Geo_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(gFieldFE,DecompFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(gFieldFE,N_gFieldFEVars,Err)
  CALL cmfe_Field_NumberOfComponentsSet(gFieldFE,Type_gFieldFEVar,N_gFieldFEComps,Err)
  DO i=1,N_gFieldFEComps
    CALL cmfe_Field_ComponentMeshComponentSet(gFieldFE,Type_gFieldFEVar,i,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(gFieldFE,U_type,1,MeshFEComp1,Err) !
  !CALL cmfe_Field_ComponentMeshComponentSet(gFieldFE,U_type,2,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(gFieldFE,U_type,3,MeshFEComp1,Err)
  ENDDO
  CALL cmfe_Field_VariableLabelSet(gFieldFE,U_type,"Geometry",Err)
  CALL cmfe_Field_CreateFinish(gFieldFE,Err)
  
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GenMesh,gFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a fibre field and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(fField,Err)
  CALL cmfe_Field_CreateStart(UN_fField,RFE,fField,Err)
  CALL cmfe_Field_TypeSet(fField,F_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(fField,DecompFE,Err)
!------------------------------------------------------------------------------------------------------------------
!   sets $ fField%GEOMETRIC_FIELD => gFieldFE                                                         |
!-----------------------------------------------------------------------------------------------------------------|
  CALL cmfe_Field_GeometricFieldSet(fField,gFieldFE,Err)                                              !
!   note: setting a geometric field for a field only works for field types: "FIELD_FIBRE_TYPE, FIELD_GENERAL_TYPE,|
!   FIELD_MATERIAL_TYPE or FIELD_GEOMETRIC_GENERAL_TYPE"                                                          |
!-----------------------------------------------------------------------------------------------------------------*
  CALL cmfe_Field_NumberOfVariablesSet(fField,N_fFieldVars,Err)
  CALL cmfe_Field_NumberOfComponentsSet(fField,Type_fFieldVar,N_fFieldComps,Err)  
  DO i=1,N_fFieldComps
   CALL cmfe_Field_ComponentMeshComponentSet(fField,Type_fFieldVar,i,MeshFEComp1,Err)  
  ENDDO
  CALL cmfe_Field_VariableLabelSet(fField,Type_fFieldVar,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(fField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a material field for Finite Elasticity and attach it to the geometric field - quadratic interpolation
  CALL cmfe_Field_Initialise(mFieldFE,Err)
  CALL cmfe_Field_CreateStart(UN_mFieldFE,RFE,mFieldFE,Err)
  CALL cmfe_Field_TypeSet(mFieldFE,Mat_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(mFieldFE,DecompFE,Err)
  CALL cmfe_Field_GeometricFieldSet(mFieldFE,gFieldFE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(mFieldFE,N_mFieldFEVars,Err)
  CALL cmfe_Field_VariableTypesSet(mFieldFE,Type_mFieldFEVar,Err)
  CALL cmfe_Field_NumberOfComponentsSet(mFieldFE,Type_mFieldFEVar(1),N_mFieldFEVar1Comps,Err)
  CALL cmfe_Field_NumberOfComponentsSet(mFieldFE,Type_mFieldFEVar(2),N_mFieldFEVar2Comps,Err)
  DO i=1,N_mFieldFEVar1Comps
    CALL cmfe_Field_ComponentInterpolationSet(mFieldFE,Type_mFieldFEVar(1),i,C_IP,Err)
  ENDDO
  CALL cmfe_Field_ComponentInterpolationSet(mFieldFE,Type_mFieldFEVar(2),1,C_IP,Err)
  CALL cmfe_Field_VariableLabelSet(mFieldFE,Type_mFieldFEVar(1),"MaterialFE",Err)
  CALL cmfe_Field_VariableLabelSet(mFieldFE,Type_mFieldFEVar(2),"Gravity",Err)
  CALL cmfe_Field_CreateFinish(mFieldFE,Err)
  !Set "scaled" Mooney-Rivlin constants c10 and c01 and others
  CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(1),SET_type,1,MAT_FE(1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(1),SET_type,2,MAT_FE(2),Err)
  CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(1),SET_type,3,MAT_FE(3),Err)
  CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(1),SET_type,4,MAT_FE(4),Err)
  CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(1),SET_type,5,0.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(1),SET_type,6,TK_lin_param,Err)
! CALL cmfe_Field_ComponentValuesInitialise(mFieldFE,Type_mFieldFEVar(2),SET_type,1,0.0_CMISSRP,Err)


! die Variablen die im FE berechnet werde nsollen (l und l punkt)
  !Create the dependent field for FE with 2 variables and * components 
  !  3-d: 3 displacement (quad interpol), 1 pressure (lin interpol) --> * = 4
  !  2-d: 2 displacement (quad interpol), 1 pressure (lin interpol) --> * = 3
  !  1-d: 1 displacement (quad interpol), 1 pressure (lin interpol) --> * = 2
  CALL cmfe_Field_Initialise(dFieldFE,Err)
  CALL cmfe_Field_CreateStart(UN_dFieldFE,RFE,dFieldFE,Err)
  CALL cmfe_Field_TypeSet(dFieldFE,GeoGen_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(dFieldFE,DecompFE,Err)
  CALL cmfe_Field_GeometricFieldSet(dFieldFE,gFieldFE,Err)
  CALL cmfe_Field_DependentTypeSet(dFieldFE,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(dFieldFE,N_dFieldFEVars,Err)
  CALL cmfe_Field_VariableTypesSet(dFieldFE,Type_dFieldFEVar,Err)
  CALL cmfe_Field_NumberOfComponentsSet(dFieldFE,Type_dFieldFEVar(1),N_dFieldFEComps,Err)
  CALL cmfe_Field_NumberOfComponentsSet(dFieldFE,Type_dFieldFEVar(2),N_dFieldFEComps,Err)
  DO i=1,N_dFieldFEComps-1
    CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(1),i,MeshFEComp1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(2),i,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(1),1,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(1),2,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(1),3,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(1),4,MeshFEComp2,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(2),1,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(2),2,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(2),3,MeshFEComp1,Err)
  !CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(2),4,MeshFEComp2,Err)
  ENDDO 
  CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(1),4,MeshFEComp2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(dFieldFE,Type_dFieldFEVar(2),4,MeshFEComp2,Err)
!  CALL cmfe_Field_ScalingTypeSet(dFieldFE,CMFE_FIELD_UNIT_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(dFieldFE,Type_dFieldFEVar(1),"DependentFE",Err)
  CALL cmfe_Field_VariableLabelSet(dFieldFE,Type_dFieldFEVar(2),"Reaction_Force",Err)
  CALL cmfe_Field_CreateFinish(dFieldFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
! die Variablen, die im FE benötigt werden von der Monodomain ("gamma")
  !Create the independent field for the active stress in the mechanics mesh
!  independent_field_auto_create = .TRUE.
  CALL cmfe_Field_Initialise(iFieldFE,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL cmfe_Field_CreateStart(UN_iFieldFE,RFE,iFieldFE,Err)
    CALL cmfe_Field_TypeSet(iFieldFE,Gen_type,Err)
    CALL cmfe_Field_MeshDecompositionSet(iFieldFE,DecompFE,Err)
    CALL cmfe_Field_GeometricFieldSet(iFieldFE,gFieldFE,Err)
    CALL cmfe_Field_DependentTypeSet(iFieldFE,CMFE_FIELD_INDEPENDENT_TYPE,Err)
    CALL cmfe_Field_NumberOfVariablesSet(iFieldFE,N_iFieldFEVars,Err)
    CALL cmfe_Field_VariableTypesSet(iFieldFE,Type_iFieldFEVar,Err)
    CALL cmfe_Field_DimensionSet(iFieldFE,Type_iFieldFEVar(1),CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(iFieldFE,Type_iFieldFEVar(1),N_iFieldFEVar1Comps,Err)
    CALL cmfe_Field_NumberOfComponentsSet(iFieldFE,Type_iFieldFEVar(2),N_iFieldFEVar2Comps,Err)
    DO i=1,N_iFieldFEVar1Comps
      CALL cmfe_Field_ComponentInterpolationSet(iFieldFE,Type_iFieldFEVar(1),i,IP_iFieldFEVar1(i),Err)
    ENDDO
    DO i=1,N_iFieldFEVar2Comps
      CALL cmfe_Field_ComponentInterpolationSet(iFieldFE,Type_iFieldFEVar(2),i,IP_iFieldFEVar2(i),Err)
    ENDDO
    CALL cmfe_Field_ComponentMeshComponentSet(iFieldFE,Type_iFieldFEVar(1),1,MeshFEComp1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(iFieldFE,Type_iFieldFEVar(1),2,MeshFEComp1,Err) !Type_iFieldFEVar(1), 3-6? Type_iFieldFEVar(2), 1-4?
    CALL cmfe_Field_DataTypeSet(iFieldFE,Type_iFieldFEVar(1),CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_DataTypeSet(iFieldFE,Type_iFieldFEVar(2),CMFE_FIELD_INTG_TYPE,Err)
    CALL cmfe_Field_VariableLabelSet(iFieldFE,Type_iFieldFEVar(1),"Active_Stress_FE",Err)
    CALL cmfe_Field_VariableLabelSet(iFieldFE,Type_iFieldFEVar(2),"subgrid_info",Err)
    CALL cmfe_Field_CreateFinish(iFieldFE,Err)
  ENDIF


  !================================================================================================================================
  !  M O N O D O M A I N
  !================================================================================================================================
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a geometric field for monodomain - quadratic interpolation
  CALL cmfe_Field_Initialise(gFieldM,Err)
  CALL cmfe_Field_CreateStart(UN_gFieldM,RM,gFieldM,Err)
  CALL cmfe_Field_TypeSet(gFieldM,Geo_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(gFieldM,DecompM,Err)
  CALL cmfe_Field_TypeSet(gFieldM,Geo_type,Err)
  CALL cmfe_Field_NumberOfVariablesSet(gFieldM,N_gFieldMVars,Err)
  CALL cmfe_Field_NumberOfComponentsSet(gFieldM,Type_gFieldMVar,N_gFieldMComps,Err)
  CALL cmfe_Field_VariableLabelSet(gFieldM,Type_gFieldMVar,"GeometryM",Err)
  CALL cmfe_Field_CreateFinish(gFieldM,Err)



  !--------------------------------------------------------------------------------------------------------------------------------
  !Create a materials field for monodomain and attach it to the geometric field - constant interpolation
  CALL cmfe_Field_Initialise(mFieldM,Err)
  CALL cmfe_Field_CreateStart(UN_mFieldM,RM,mFieldM,Err)
  CALL cmfe_Field_TypeSet(mFieldM,Mat_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(mFieldM,DecompM,Err)
  CALL cmfe_Field_GeometricFieldSet(mFieldM,gFieldM,Err)
  CALL cmfe_Field_NumberOfVariablesSet(mFieldM,N_mFieldMVars,Err)
  CALL cmfe_Field_NumberOfComponentsSet(mFieldM,Type_mFieldMVar,N_mFieldMComps,Err)
  DO i=1,N_mFieldMComps
    CALL cmfe_Field_ComponentInterpolationSet(mFieldM,Type_mFieldMVar,i,IP_mFieldMVar(i),Err)
  ENDDO
  CALL cmfe_Field_VariableLabelSet(mFieldM,Type_mFieldMVar,"MaterialM",Err)
  CALL cmfe_Field_CreateFinish(mFieldM,Err)
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(mFieldM,Type_mFieldMVar,SET_type,1,Am,Err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(mFieldM,Type_mFieldMVar,SET_type,2,Cm_fast,Err)
  !Set Conductivity
  CALL cmfe_Field_ComponentValuesInitialise(mFieldM,Type_mFieldMVar,SET_type,3,CONDUCTIVITY,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the dependent field for monodomain with 2 variables and 1 components 
  CALL cmfe_Field_Initialise(dFieldM,Err)
  CALL cmfe_Field_CreateStart(UN_dFieldM,RM,dFieldM,Err)
  CALL cmfe_Field_TypeSet(dFieldM,Gen_type,Err)
  CALL cmfe_Field_MeshDecompositionSet(dFieldM,DecompM,Err)
  CALL cmfe_Field_GeometricFieldSet(dFieldM,gFieldM,Err)
  CALL cmfe_Field_DependentTypeSet(dFieldM,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(dFieldM,N_dFieldMVars,Err)
  CALL cmfe_Field_VariableTypesSet(dFieldM,Type_dFieldMVar,Err)
  CALL cmfe_Field_NumberOfComponentsSet(dFieldM,Type_dFieldMVar(1),N_dFieldMVar1Comps,Err)
  CALL cmfe_Field_NumberOfComponentsSet(dFieldM,Type_dFieldMVar(2),N_dFieldMVar2Comps,Err)
  !additionally the V_Var_Type (Type_dFieldMVar(3)) with 3 components for the 3-d position of the geometry:
  CALL cmfe_Field_NumberOfComponentsSet(dFieldM,Type_dFieldMVar(3),N_dFieldMVar3Comps,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(dFieldM,Type_dFieldMVar(1),1,MeshMComp1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(dFieldM,Type_dFieldMVar(2),1,MeshMComp1,Err)
!  CALL cmfe_Field_ComponentMeshComponentSet(dFieldM,Type_dFieldMVar(3),1,MeshMComp1,Err)
  CALL cmfe_Field_DimensionSet(dFieldM,Type_dFieldMVar(1),CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_DimensionSet(dFieldM,Type_dFieldMVar(2),CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(dFieldM,Type_dFieldMVar(1),"Vm",Err)
  CALL cmfe_Field_VariableLabelSet(dFieldM,Type_dFieldMVar(2),"dVm/dt",Err)
  CALL cmfe_Field_VariableLabelSet(dFieldM,Type_dFieldMVar(3),"GeometryM3D",Err)
  CALL cmfe_Field_CreateFinish(dFieldM,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the independent field for the active stress in the electrics mesh
!  independent_field_auto_create = .TRUE.
  CALL cmfe_Field_Initialise(iFieldM,Err)
  IF(.NOT. independent_field_auto_create) THEN
    CALL cmfe_Field_CreateStart(UN_iFieldM,RM,iFieldM,Err)
    CALL cmfe_Field_TypeSet(iFieldM,Gen_type,Err)
    CALL cmfe_Field_MeshDecompositionSet(iFieldM,DecompM,Err)
    CALL cmfe_Field_GeometricFieldSet(iFieldM,gFieldM,Err)
    CALL cmfe_Field_DependentTypeSet(iFieldM,CMFE_FIELD_INDEPENDENT_TYPE,Err)
    CALL cmfe_Field_NumberOfVariablesSet(iFieldM,N_iFieldMVars,Err)
    CALL cmfe_Field_VariableTypesSet(iFieldM,Type_iFieldMVar,Err)
    !first variable:   U_type
    CALL cmfe_Field_DataTypeSet(iFieldM,Type_iFieldMVar(1),CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_DimensionSet(iFieldM,Type_iFieldMVar(1),CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(iFieldM,Type_iFieldMVar(1),N_iFieldMVar1Comps,Err)
    DO i=1,N_iFieldMVar1Comps
      CALL cmfe_Field_ComponentInterpolationSet(iFieldM,Type_iFieldMVar(1),i,IP_iFieldMVar1(i),Err)
    ENDDO
    CALL cmfe_Field_VariableLabelSet(iFieldM,Type_iFieldMVar(1),"Active_Stress_M",Err)
    !second variable:   V_type -- 1) motor unit number   2) fibre type   3) fibre number   4) nearest Gauss point   5) in element number (LOCAL NODE NUMBERING!!!)
    CALL cmfe_Field_DataTypeSet(iFieldM,Type_iFieldMVar(2),CMFE_FIELD_INTG_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(iFieldM,Type_iFieldMVar(2),N_iFieldMVar2Comps,Err)
    DO i=1,N_iFieldMVar2Comps
      CALL cmfe_Field_ComponentInterpolationSet(iFieldM,Type_iFieldMVar(2),i,IP_iFieldMVar2(i),Err)
    ENDDO
    CALL cmfe_Field_VariableLabelSet(iFieldM,Type_iFieldMVar(2),"fibre_info",Err)
    !third variable:   FIELD_U1_VARIABLE_TYPE -- 1) sarcomere half length   2) inital sarcomere half length   3) initial node distance   4) sarcomere half length when activation starts
    CALL cmfe_Field_DataTypeSet(iFieldM,Type_iFieldMVar(3),CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(iFieldM,Type_iFieldMVar(3),N_iFieldMVar3Comps,Err)
    DO i=1,N_iFieldMVar3Comps 
      CALL cmfe_Field_ComponentInterpolationSet(iFieldM,Type_iFieldMVar(3),i,IP_iFieldMVar3(i),Err)
    ENDDO
    CALL cmfe_Field_VariableLabelSet(iFieldM,Type_iFieldMVar(3),"sarcomere_half_length",Err)
    !fourth variable:   FIELD_U2_VARIABLE_TYPE -- 1) old node distance   2) maximum contraction velocity   3) relative contraction velocity
    CALL cmfe_Field_DataTypeSet(iFieldM,Type_iFieldMVar(4),CMFE_FIELD_DP_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(iFieldM,Type_iFieldMVar(4),N_iFieldMVar4Comps,Err)
    DO i=1,N_iFieldMVar4Comps
      CALL cmfe_Field_ComponentInterpolationSet(iFieldM,Type_iFieldMVar(4),i,IP_iFieldMVar4(i),Err)
    ENDDO
    CALL cmfe_Field_VariableLabelSet(iFieldM,Type_iFieldMVar(4),"contraction_velocity",Err)
    CALL cmfe_Field_CreateFinish(iFieldM,Err)
    !todo : pack labels into vector, 'automise' the field initialisations
  ENDIF
!=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?||

  
  !================================================================================================================================
  !  EQUATIONS SET

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for Finite Elasticity
  CALL cmfe_Field_Initialise(esFieldFE,Err)
  CALL cmfe_EquationsSet_Initialise(EqusSetFE,Err)
  CALL cmfe_EquationsSet_CreateStart(UN_EquSetFE,RFE,fField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
   & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE], & 
   & UN_esFieldFE,esFieldFE,EqusSetFE,Err)
  CALL cmfe_EquationsSet_CreateFinish(EqusSetFE,Err)

  !Create the equations set dependent field variables for Finite Elasticity
  CALL cmfe_EquationsSet_DependentCreateStart(EqusSetFE,UN_dFieldFE,dFieldFE,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EqusSetFE,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EqusSetFE,UN_iFieldFE,iFieldFE,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EqusSetFE,Err)

  !Create the equations set materials field variables for Finite Elasticity
  CALL cmfe_EquationsSet_MaterialsCreateStart(EqusSetFE,UN_mFieldFE,mFieldFE,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EqusSetFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations_set for monodomain
  CALL cmfe_Field_Initialise(esFieldM,Err)
  CALL cmfe_EquationsSet_Initialise(EqusSetM,Err)
  !Set the equations set to be a Monodomain equations set
  !> \todo solve the monodomain problem on the fibre field rather than on the geometric field: GeometricField <--> fField
  CALL cmfe_EquationsSet_CreateStart(UN_EquSetM,RM,gFieldM,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
   & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE], &
   & UN_esFieldM,esFieldM,EqusSetM,Err)
  CALL cmfe_EquationsSet_CreateFinish(EqusSetM,Err)

  !Create the equations set dependent field variables for monodomain
  CALL cmfe_EquationsSet_DependentCreateStart(EqusSetM,UN_dFieldM,dFieldM,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EqusSetM,Err)

  !Create the equations set independent field variable for the active stress component for Finite Elasticity
  CALL cmfe_EquationsSet_IndependentCreateStart(EqusSetM,UN_iFieldM,iFieldM,Err)
  CALL cmfe_EquationsSet_IndependentCreateFinish(EqusSetM,Err)

  !Create the equations set materials field variables for monodomain
  CALL cmfe_EquationsSet_MaterialsCreateStart(EqusSetM,UN_mFieldM,mFieldM,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EqusSetM,Err)

 
  !--------------------------------------------------------------------------------------------------------------------------------
  !UPDATE THE INDEPENDENT FIELD iFieldM
  !first variable
  !  components: 
  !    1) active stress
!  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U_type,SET_type,2,0.0_CMISSRP,Err)
!  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U_type,SET_type,3,0.0_CMISSRP,Err)
!  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U_type,SET_type,4,0.0_CMISSRP,Err)
  !
  !second variable
  !  components: 
  !    1) motor unit number
  !    2) fibre type
  !    3) fibre number
  !    4) nearest Gauss point
  !    5) in element number (LOCAL NODE NUMBERING!!!)
  !
  !init the motor unit number and fibre type to 1
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,V_type,SET_type,1,1,Err) !mu_nr=1
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,V_type,SET_type,2,1,Err) !Ftype=1
  !init the fibre number, the nearest Gauss point info and the inElem info to 0
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,V_type,SET_type,3,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,V_type,SET_type,4,0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,V_type,SET_type,5,0,Err) !(LOCAL NODE NUMBERING!!!)
  !third variable:
  !  components:
  !    1) sarcomere half length
  !    2) initial sarcomere half length
  !    3) initial node distance
  !    4) sarcomere half length when activation starts
! put upwards  stretch_sarcolength_ratio=1.0_CMISSRP/1.0_CMISSRP
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U1_type,SET_type,2, &
!   & 1.0_CMISSRP,Err)
   & stretch_sarcolength_ratio,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U1_type,SET_type,3, &
   & LENGTH/N_FibresInSeries/(N_NodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U1_type,SET_type,4, &
   & 1.0_CMISSRP,Err)
  !fourth variable:
  !  components:
  !    1) old node distance
  !    2) maximum contraction velocity
  !    3) relative contraction velocity
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,1, &
   & LENGTH/N_FibresInSeries/(N_NodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,2, &
   & Vmax/N_FibresInSeries/(N_NodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,3, &
   & 0.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,4, &
   & LENGTH/N_FibresInSeries/(N_NodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,5, &
   & LENGTH/N_FibresInSeries/(N_NodesPerFibre-1),Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,6, &
   & LENGTH/N_FibresInSeries/(N_NodesPerFibre-1),Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !UPDATE THE INDEPENDENT FIELD iFieldFE
  !second variable of iFieldFE
  !  components:
  !    1) number of nodes in Xi(1) direction per element
  !    2) number of nodes in Xi(2) direction per element
  !    3) number of nodes in Xi(3) direction per element
  !    4) beginning of fibres in this FE element??? 1=yes, 0=no
  !
  !initialise as if the fibres would not start in any element, and adjust below
  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,1,N_NodesInXi1-1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,2,N_NodesInXi2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,3,N_NodesInXi3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,4,0,Err)

  !fibres are starting in elements 1,4,7,10,...
  DO elem_idx=1,N_ElementsFE,N_GlXElements/N_FibresInSeries
    CALL cmfe_Decomposition_ElementDomainGet(DecompFE,elem_idx,ElementDomain,Err)
    IF(ElementDomain==k_CNode) THEN
      !fibres begin in this element
      CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type,elem_idx,4,1,Err)
      CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type,elem_idx,1,&
           & N_NodesInXi1,Err)
    ENDIF
  ENDDO

!!  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,4,0,Err)
!  !fibres start in all elements
!  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,4,1,Err)
!  !adjust numer of nodes in each element in Xi1 direction
!  CALL cmfe_Field_ComponentValuesInitialise(iFieldFE,V_type,SET_type,1, &
!   & N_NodesInXi1,Err)

!  !fibres are starting only in element 1
!  elem_idx=1
!  CALL cmfe_Decomposition_ElementDomainGet(DecompFE,elem_idx,ElementDomain,Err)
!  IF(ElementDomain==k_CNode) THEN
!    !fibres begin in this element
!    CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type,elem_idx,4,1,Err) 
!    CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type, &
!     & elem_idx,1,N_NodesInXi1,Err)
!  ENDIF

!  !fibres are starting in elements 1,3,5, and 7
!  DO elem_idx=1,N_ElementsFE,2
!    CALL cmfe_Decomposition_ElementDomainGet(DecompFE,elem_idx,ElementDomain,Err)
!    IF(ElementDomain==k_CNode) THEN
!      !fibres begin in this element
!      CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type,elem_idx,4,1,Err) 
!      CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type, &
!       & elem_idx,1,N_NodesInXi1,Err)
!    ENDIF
!  ENDDO

!  !fibres are not starting in elements 2,4,6, and 8
!  DO elem_idx=2,N_ElementsFE,2
!    CALL cmfe_Decomposition_ElementDomainGet(DecompFE,elem_idx,ElementDomain,Err)
!    IF(ElementDomain==k_CNode) THEN
!      CALL cmfe_Field_ParameterSetUpdateElement(iFieldFE,V_type,SET_type, &
!       & elem_idx,4,0,Err) !fibres do not begin in this element
!    ENDIF
!  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the equations set equations for monodomain
  CALL cmfe_Equations_Initialise(EqusM,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EqusSetM,EqusM,Err)
  CALL cmfe_Equations_SparsityTypeSet(EqusM,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EqusM,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EqusSetM,Err)

  !Create the equations set equations for Finite Elasticity
  CALL cmfe_Equations_Initialise(EqusFE,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EqusSetFE,EqusFE,Err)
  CALL cmfe_Equations_SparsityTypeSet(EqusFE,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EqusFE,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EqusSetFE,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML environment
  CALL cmfe_CellML_Initialise(Cml,Err)
  CALL cmfe_CellML_CreateStart(UN_Cml,RM,Cml,Err)
  !Import the Shorten et al. 2007 model from a file
  CALL cmfe_CellML_ModelImport(Cml,filename,idx_CmlSH,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !,- to set from this side
  CALL cmfe_CellML_VariableSetAsKnown(Cml,idx_CmlSH,"Aliev_Panfilov/I_HH",Err)
  CALL cmfe_CellML_VariableSetAsKnown(Cml,idx_CmlSH,"Razumova/l_hs",Err)
  CALL cmfe_CellML_VariableSetAsKnown(Cml,idx_CmlSH,"Razumova/rel_velo",Err)
!
!  CALL cmfe_CellML_VariableSetAsKnown(Cml,shortenModelIndex2,"wal_environment/I_HH",Err)
!  CALL cmfe_CellML_VariableSetAsKnown(Cml,shortenModelIndex2,"razumova/L_S",Err)
!  CALL cmfe_CellML_VariableSetAsKnown(Cml,shortenModelIndex2,"razumova/rel_velo",Err)
  !,- to get from the CellML side
!  CALL cmfe_CellML_VariableSetAsWanted(Cml,idx_CmlSH,"wal_environment/I_T",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(Cml,idx_CmlSH,"wal_environment/I_ionic_s",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(Cml,idx_CmlSH,"wal_environment/I_ionic_t",Err)
  !
  !NOTE: If an INTERMEDIATE (or ALGEBRAIC) variable should be used in a mapping, it has to be set as known or wanted first!
  !,  --> set "razumova/stress" as wanted!
  !,  --> no need to set "wal_environment/vS" since all STATE variables are automatically set as wanted! 
  CALL cmfe_CellML_VariableSetAsWanted(Cml,idx_CmlSH,"Razumova/ActiveStress",Err)
  CALL cmfe_CellML_VariableSetAsWanted(Cml,idx_CmlSH,"Razumova/Activation",Err)
!  CALL cmfe_CellML_VariableSetAsWanted(Cml,shortenModelIndex2,"razumova/stress",Err)
  !,- and override constant parameters without needing to set up fields
  !> \todo Need to allow parameter values to be overridden for the case when user has non-spatially varying parameter value.
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(Cml,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(Cml,Err)
  !Map the sarcomere half length L_S
  CALL cmfe_CellML_CreateFieldToCellMLMap(Cml,iFieldM,U1_type,1,SET_type, &
   & idx_CmlSH,"Razumova/l_hs",SET_type,Err)
  !Map the sarcomere relative contraction velocity
  CALL cmfe_CellML_CreateFieldToCellMLMap(Cml,iFieldM,U2_type,3,SET_type, &
   & idx_CmlSH,"Razumova/rel_velo",SET_type,Err)
  !Map the transmembrane voltage V_m
  CALL cmfe_CellML_CreateFieldToCellMLMap(Cml,dFieldM,U_type,1,SET_type, &
    & idx_CmlSH,"Aliev_Panfilov/V_m",SET_type,Err)
!   & idx_CmlSH,"wal_environment/vS",SET_type,Err)
!  CALL cmfe_CellML_CreateCellMLToFieldMap(Cml,idx_CmlSH,"wal_environment/vS",SET_type, &
  CALL cmfe_CellML_CreateCellMLToFieldMap(Cml,idx_CmlSH,"Aliev_Panfilov/V_m",SET_type, &
   & dFieldM,U_type,1,SET_type,Err)
  !Map the active stress
  CALL cmfe_CellML_CreateCellMLToFieldMap(Cml,idx_CmlSH,"Razumova/ActiveStress",SET_type, &
   & iFieldM,U_type,1,SET_type,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(Cml,idx_CmlSH,"Razumova/Activation",SET_type, &
   & iFieldM,U_type,6,SET_type,Err)

!  CALL cmfe_CellML_CreateFieldToCellMLMap(Cml,iFieldM,U1_type,1,SET_type, &
!   & shortenModelIndex2,"razumova/L_S",SET_type,Err)
!  !Map the sarcomere relative contraction velocity
!  CALL cmfe_CellML_CreateFieldToCellMLMap(Cml,iFieldM,U2_type,3,SET_type, &
!   & shortenModelIndex2,"razumova/rel_velo",SET_type,Err)
!  !Map the transmembrane voltage V_m
!  CALL cmfe_CellML_CreateFieldToCellMLMap(Cml,dFieldM,U_type,1,SET_type, &
!   & shortenModelIndex2,"wal_environment/vS",SET_type,Err)
!  CALL cmfe_CellML_CreateCellMLToFieldMap(Cml,shortenModelIndex2,"wal_environment/vS",SET_type, &
!   & dFieldM,U_type,1,SET_type,Err)
!  !Map the active stress
!  CALL cmfe_CellML_CreateCellMLToFieldMap(Cml,shortenModelIndex2,"razumova/stress",SET_type, &
!   & iFieldM,U_type,1,SET_type,Err)

  CALL cmfe_CellML_FieldMapsCreateFinish(Cml,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Initialise dependent field for monodomain
  !> \todo - get V_m initialial value.
  CALL cmfe_Field_ComponentValuesInitialise(dFieldM,U_type,SET_type,1, &
   !& -79.974_CMISSRP,Err)
   & -0.0_CMISSRP,Err)
  
  !Initialise dependent field for Finite Elasticity from undeformed geometry and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(gFieldFE,U_type, &
   & SET_type,1,dFieldFE,U_type,SET_type,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(gFieldFE,U_type, & 
   & SET_type,2,dFieldFE,U_type,SET_type,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(gFieldFE,U_type, &
   & SET_type,3,dFieldFE,U_type,SET_type,3,Err)
! put upwards  INIT_PRESSURE=-2.0_CMISSRP*MAT_FE(2)-MAT_FE(1)
  CALL cmfe_Field_ComponentValuesInitialise(dFieldFE,U_type,SET_type,4,INIT_PRESSURE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the CellML models field
  CALL cmfe_Field_Initialise(msFieldCml,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(Cml,UN_msFieldCml,msFieldCml,Err)
  CALL cmfe_CellML_ModelsFieldCreateFinish(Cml,Err)

  !Set up the models field
  CALL cmfe_Field_ComponentValuesInitialise(msFieldCml,U_type,SET_type,1,idx_CmlSH,Err)

!  DO k_Node=N_NodesPerFibre/2,N_NodesM,N_NodesPerFibre
!    CALL cmfe_Decomposition_NodeDomainGet(DecompM,k_Node,1,NodeDomain,Err)
!    IF(NodeDomain==k_CNode) THEN
!      CALL cmfe_Field_ParameterSetUpdateNode(msFieldCml,U_type, &
!        & SET_type,1,1,k_Node,1,shortenModelIndex2,Err)
!    ENDIF
!  ENDDO

!  CALL cmfe_Field_ParameterSetUpdateStart(msFieldCml,U_type,SET_type,Err)
!  CALL cmfe_Field_ParameterSetUpdateFinish(msFieldCml,U_type,SET_type,Err)

  !Create the CellML state field
  CALL cmfe_Field_Initialise(sFieldCml,Err)
  CALL cmfe_CellML_StateFieldCreateStart(Cml,UN_sFieldCml,sFieldCml,Err)
  CALL cmfe_CellML_StateFieldCreateFinish(Cml,Err)

  !Create the CellML intermediate field
  CALL cmfe_Field_Initialise(imFieldCml,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(Cml,UN_imFieldCml,imFieldCml,Err)
  CALL cmfe_CellML_IntermediateFieldCreateFinish(Cml,Err)
  
  !Create the CellML parameters field
  CALL cmfe_Field_Initialise(pFieldCml,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(Cml,UN_pFieldCml,pFieldCml,Err)
  CALL cmfe_CellML_ParametersFieldCreateFinish(Cml,Err)
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(UN_Problem,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
   & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE],Problem,Err)
!  CALL cmfe_Problem_SpecificationSet(Problem,CMFE_PROBLEM_MULTI_PHYSICS_CLASS,CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE, &
!   & CMFE_PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)

  !set the main control loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopMain,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoopMain,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopMain,'MAIN_TIME_LOOP',Err)
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,0.0_CMISSRP,ELASTICITY_TIME_STEP,ELASTICITY_TIME_STEP,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,OUTPUT_FREQUENCY,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopMain,CMFE_CONTROL_LOOP_TIMING_OUTPUT,Err) !DO NOT CHANGE!!!

  !set the monodomain loop (time loop type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[k_ControlLoopM,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopM,'MONODOMAIN_TIME_LOOP',Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoopM,0.0_CMISSRP,ELASTICITY_TIME_STEP,PDE_TIME_STEP,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoopM,CMFE_CONTROL_LOOP_NO_OUTPUT,Err)

  !set the finite elasticity loop (simple type)
  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[k_ControlLoopFE,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_TypeSet(ControlLoopFE,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err) ! tomo
!  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoopFE,'ELASTICITY_LOOP',Err)

  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[k_ControlLoopM,CMFE_CONTROL_LOOP_NODE], &
   & idx_SolDAE,SolverDAE,Err)
  CALL cmfe_Solver_DAETimeStepSet(SolverDAE,ODE_TIME_STEP,Err)
  !> \todo - solve the CellML equations on the GPU for efficiency (later)
  !CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EXTERNAL,Err) 
  CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Create the parabolic solver
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_Problem_SolverGet(Problem,[k_ControlLoopM,CMFE_CONTROL_LOOP_NODE], &
   & idx_SolPar,SolverParabolic,Err)
  CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Create the Finte Elasticity solver
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_Solver_Initialise(LinearSolverFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[k_ControlLoopFE,CMFE_CONTROL_LOOP_NODE], &
   & idx_SolFE,SolverFE,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(SolverFE,500,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(SolverFE,1.E-6_CMISSRP,Err) !6
  CALL cmfe_Solver_NewtonSolutionToleranceSet(SolverFE,2.E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
!  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)

  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)

  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_CellMLEquations_Initialise(EqusCml,Err)
  CALL cmfe_Problem_SolverGet(Problem,[k_ControlLoopM,CMFE_CONTROL_LOOP_NODE],idx_SolDAE,SolverDAE,Err)
  CALL cmfe_Solver_CellMLEquationsGet(SolverDAE,EqusCml,Err)
  ! adds cellML environment to cellMLequations and sets index of environment:
  CALL cmfe_CellMLEquations_CellMLAdd(EqusCml,Cml,idx_Cml,Err)

  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_SolverEquations_Initialise(SolEqusM,Err)
  CALL cmfe_Problem_SolverGet(Problem,[k_ControlLoopM,CMFE_CONTROL_LOOP_NODE], &
   & idx_SolPar,SolverParabolic,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverParabolic,SolEqusM,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolEqusM,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolEqusM,CMFE_SOLVER_FULL_MATRICES,Err)  
  CALL cmfe_SolverEquations_EquationsSetAdd(SolEqusM,EqusSetM,indx_EqusSetM,Err)

  !Create the problem solver Finite Elasticity equations
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_SolverEquations_Initialise(SolEqusFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[k_ControlLoopFE,CMFE_CONTROL_LOOP_NODE],idx_SolFE,SolverFE,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverFE,SolEqusFE,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolEqusFE,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolEqusFE,CMFE_SOLVER_FULL_MATRICES,Err)
  !Arguments:  solverEqus - The ~ to add the equs set for.  equsSet -The ~ to add.   equsSetIndex -On return, the index of the added equsset in the solver equs.
  CALL cmfe_SolverEquations_EquationsSetAdd(SolEqusFE,EqusSetFE,idx_EqusSetFE,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !boundary conditions

  !Prescribe boundary conditions for monodomain
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolEqusM,BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolEqusM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolEqusFE,BoundaryConditionsFE,Err)

  CALL cmfe_GeneratedMesh_SurfaceGet(GenMesh,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GenMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GenMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GenMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 100% x displacement
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    k_Node=LeftSurfaceNodes(node_idx)
    ! #'1' is the meshComponentNumber - The user number of the mesh component to get the domain for. #'NodeDomain' - On return, the computational domain of the node.
    CALL cmfe_Decomposition_NodeDomainGet(DecompFE,k_Node,1,NodeDomain,Err)
    IF(NodeDomain==k_CNode) THEN ! (to see whether the node is to be treated by 'this' CNode)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,dFieldFE,U_type,1,1,k_Node,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    k_Node=RightSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompFE,k_Node,1,NodeDomain,Err)
    IF(NodeDomain==k_CNode) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,dFieldFE,U_type,1,1,k_Node,1, &
!        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,LENGTH*1.0_CMISSRP,Err) ! To change the initial half-sarcomere length. To create cases similar to Rode. 
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,LENGTH*1.0_CMISSRP/stretch_sarcolength_ratio,Err) ! To change the initial half-sarcomere length. To create cases similar to Rode. 
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    k_Node=FrontSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompFE,k_Node,1,NodeDomain,Err)
    IF(NodeDomain==k_CNode) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,dFieldFE,U_type,1,1,k_Node,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    k_Node=BottomSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompFE,k_Node,1,NodeDomain,Err)
    IF(NodeDomain==k_CNode) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,dFieldFE,U_type,1,1,k_Node,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO


  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolEqusFE,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Calculate the bioelectrics geometric field 
  CALL cmfe_ControlLoop_Initialise(ControlLoopM,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[k_ControlLoopM,CMFE_CONTROL_LOOP_NODE],ControlLoopM,Err)
  CALL cmfe_BioelectricsFiniteElasticity_UpdateGeometricField(ControlLoopM,.TRUE.,Err)

  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,3, &
   & 0.0_CMISSRP,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RM,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"titinExample_M","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"titinExample_M","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RFE,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"titinExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"titinExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF



  !--------------------------------------------------------------------------------------------------------------------------------
  !Solve the problem


  !Solve the problem -- bring to new length before applying the stimulus
  WRITE(*,'(A)') "Start solve before stimulation"
  CALL cmfe_Problem_Solve(Problem,Err)

  !reset the relative contraction velocity to 0
  CALL cmfe_Field_ComponentValuesInitialise(iFieldM,U2_type,SET_type,3, &
   & 0.0_CMISSRP,Err)

  CALL cmfe_ControlLoop_Initialise(ControlLoopFE,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[k_ControlLoopFE,CMFE_CONTROL_LOOP_NODE],ControlLoopFE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoopFE,1,Err)

  CALL cmfe_Field_ParameterSetUpdateConstant(mFieldFE,U_type,SET_type,5,P_max,Err)


! no change for BCs -- fix at this length!!!




  !Set the Stimulus for monodomain at the middle of the fibres
  CALL cmfe_CellML_FieldComponentGet(Cml,idx_CmlSH,CMFE_CELLML_PARAMETERS_FIELD, &
   & "Aliev_Panfilov/I_HH",stimcomponent,Err) !"wal_environment/I_HH",stimcomponent,Err)

  !update the sarcomere stretch at activation
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(iFieldM,U1_type, &
   & SET_type,1,iFieldM,U1_type,SET_type,4,Err)



  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
	STIM_VALUE=90.0_CMISSRP !setzen oder nicht? wie wird das cellml model richtig (überhaupt) getriggert?

	k_Node=(N_NodesPerFibre+1)/2
  DO WHILE(k_Node<N_NodesM)
    CALL cmfe_Decomposition_NodeDomainGet(DecompM,k_Node,1,NodeDomain,Err)
    IF(NodeDomain==k_CNode) CALL cmfe_Field_ParameterSetUpdateNode(pFieldCml, &
     & U_type,SET_type,1,1,k_Node,stimcomponent,STIM_VALUE,Err)
    k_Node=k_Node+N_NodesPerFibre
  ENDDO
  
 ! variable 'time' was never set until now!
  CALL cmfe_ControlLoop_TimesSet(ControlLoopMain,time,time+STIM_STOP,ELASTICITY_TIME_STEP,Err)

  !Solve the problem for the stimulation time
  WRITE(*,'(A)') "Start solve with stimulation"
  CALL cmfe_Problem_Solve(Problem,Err)


  
  !--------------------------------------------------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------------------------------------------------


  
  !--------------------------------------------------------------------------------------------------------------------------------
!  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RM,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"titinExample_M","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"titinExample_M","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(RFE,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"titinExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"titinExample_FE","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF
  
  !--------------------------------------------------------------------------------------------------------------------------------
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  if(N_CNodes==1) then
    open(unit=888,file="times.out",iostat=stat)
    write(888,*) 'single processor'
!!! time elapsed
    total = etime(elapsed)
    write(888,*) 'End: total=', total, ' user=', elapsed(1), ' system=', elapsed(2)
!!!
    close(unit=888)
  else
    open(unit=888,file="times.out",iostat=stat,access='append')
    write(888,*) 'two processors'
!!! time elapsed
    total = etime(elapsed)
    write(888,*) 'End: total=', total, ' user=', elapsed(1), ' system=', elapsed(2)
!!!
    close(unit=888)
  endif
  
  STOP
  
END PROGRAM TITINEXAMPLE
