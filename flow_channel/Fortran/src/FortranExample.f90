!> \file
!> \author Andreas Hessenthaler
!> \brief This is an example program to solve a Navier-Stokes equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Andreas Hessenthaler
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>
!> Main program
PROGRAM NavierStokesExample

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


  !Test program parameters
  REAL(CMISSRP), PARAMETER :: StartTime=0.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: StopTime=0.1_CMISSRP
  REAL(CMISSRP), PARAMETER :: TimeIncrement=0.001_CMISSRP
  REAL(CMISSRP), PARAMETER :: Width=10.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: Length=2.0_CMISSRP
  REAL(CMISSRP) :: Mu=1.0_CMISSRP
  REAL(CMISSRP) :: Rho=1.0_CMISSRP
  REAL(CMISSRP) :: vx,vy,x,y
  REAL(CMISSRP) :: vxmax=1.0_CMISSRP
  INTEGER(CMISSIntg) :: VelocityInterpolationType
  INTEGER(CMISSIntg) :: PressureInterpolationType
  INTEGER(CMISSIntg) :: PressureMeshComponent
  INTEGER(CMISSIntg) :: NumberOfGaussXi
  INTEGER(CMISSIntg) :: ScalingType
  INTEGER(CMISSIntg), PARAMETER :: NumberOfLoadIncrements=2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: VelocityBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldSourceUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program variables
  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,component_idx,deriv_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:), RightSurfaceNodes(:), TopSurfaceNodes(:), BottomSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi, RightNormalXi, TopNormalXi, BottomNormalXi
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument

  !CMISS variables
  TYPE(cmfe_BasisType) :: VelocityBasis,PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,SourceField,EquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,NonlinearSolver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  CALL cmfe_OutputSetOn("Channel",Err)

  !Read in arguments and overwrite default values
  !Usage: NavierStokesExample [Velocity Interpolation Type] [X elements] [Y elements] [Z elements] [Scaling Type]
  !Defaults:
  VelocityInterpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  NumberGlobalXElements=3
  NumberGlobalYElements=2
  ScalingType=CMFE_FIELD_NO_SCALING

  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalXElements
    IF(NumberGlobalXElements<1) CALL HANDLE_ERROR("Invalid number of X elements.")
  ENDIF
  IF(NumberOfArguments >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalYElements
    IF(NumberGlobalYElements<1) CALL HANDLE_ERROR("Invalid number of Y elements.")
  ENDIF
  IF(NumberOfArguments >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(CommandArgument(1:ArgumentLength),*) Mu
    IF(NumberGlobalYElements<0) CALL HANDLE_ERROR("Invalid parameter mu.")
  ENDIF
  IF(NumberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(CommandArgument(1:ArgumentLength),*) Rho
    IF(NumberGlobalYElements<0) CALL HANDLE_ERROR("Invalid parameter rho.")
  ENDIF
  IF(NumberOfArguments >= 5) THEN
    CALL GET_COMMAND_ARGUMENT(5,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 5.")
    READ(CommandArgument(1:ArgumentLength),*) vxmax
    IF(NumberGlobalYElements<0) CALL HANDLE_ERROR("Invalid parameter vxmax.")
  ENDIF
  NumberOfGaussXi=3
  PressureMeshComponent=2
  WRITE(*,'("Interpolation: ", 3 i3)') VelocityInterpolationType, PressureInterpolationType
  WRITE(*,'("Elements: ", 2 i3)') NumberGlobalXElements,NumberGlobalYElements

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis function for Velocity
  CALL cmfe_Basis_Initialise(VelocityBasis,Err)
  CALL cmfe_Basis_CreateStart(VelocityBasisUserNumber,VelocityBasis,Err)
  CALL cmfe_Basis_TypeSet(VelocityBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(VelocityBasis,2,Err)
  CALL cmfe_Basis_InterpolationXiSet(VelocityBasis,[VelocityInterpolationType,VelocityInterpolationType],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(VelocityBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL cmfe_Basis_CreateFinish(VelocityBasis,Err)

  !Basis for pressure
  CALL cmfe_Basis_Initialise(PressureBasis,Err)
  CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(PressureBasis,2,Err)
  CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
  CALL cmfe_Basis_CreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[VelocityBasis,PressureBasis],Err)
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Length,Width],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    &  CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    &  CMFE_EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  DO component_idx=1,2
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3, &
    & PressureMeshComponent,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3, &
    & PressureMeshComponent,Err)
!  CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3, &
!    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!  CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3, &
!    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ScalingTypeSet(DependentField,ScalingType,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"MuRho",Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set viscosity mu and density rho
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,Mu,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,Rho,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field as zero
  CALL cmfe_Field_ComponentValuesInitialise(DependentField, &
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField, &
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField, &
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 3,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber, &
    & [CMFE_PROBLEM_FLUID_MECHANICS_CLASS, &
    &  CMFE_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    &  CMFE_PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE], &
    & Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,StartTime,StopTime,TimeIncrement,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,1,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_DynamicThetaSet(Solver,1.0_CMISSRP,Err)
  CALL cmfe_Solver_DynamicNonlinearSolverGet(Solver,NonlinearSolver,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolver, &
    & CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolver,CMFE_SOLVER_NO_OUTPUT,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  CALL cmfe_Solver_OutputTypeSet(LinearSolver,CMFE_SOLVER_NO_OUTPUT,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Solver_LibraryTypeSet(LinearSolver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BACK_SURFACE,BottomSurfaceNodes,BottomNormalXi,Err)

  !Flow at x=0
  DO node_idx=1,SIZE(TopSurfaceNodes,1)
    NodeNumber=TopSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode( &
        & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
        & 1,1,NodeNumber,1,x,Err)
      CALL cmfe_Field_ParameterSetGetNode( &
        & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
        & 1,1,NodeNumber,2,y,Err)
      vy = -4.0_CMISSRP*x*(x-Length)/Length/Length*vxmax
      vx = 0.0_CMISSRP
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 1,CMFE_BOUNDARY_CONDITION_FIXED,vx,Err)
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 2,CMFE_BOUNDARY_CONDITION_FIXED,vy,Err)
    ENDIF
  ENDDO
  !Flow at y=0
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Flow at y=width
  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Outflow at x=length
  DO node_idx=1,SIZE(BottomSurfaceNodes,1),2
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"Channel","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"Channel","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HANDLE_ERROR

END PROGRAM NavierStokesExample

