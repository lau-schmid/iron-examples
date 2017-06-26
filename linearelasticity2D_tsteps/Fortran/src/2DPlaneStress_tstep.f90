! \file
!> \author Chris Bradley
!> \brief This is an example program to solve a linear elasticity equation using OpenCMISS calls.
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
!> Contributor(s):
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

!> \example LinearElasticity/src/LinearElasticityExample.f90
!! Example program to solve a linear elasticity equation using OpenCMISS calls.
!<

!> Main program
PROGRAM LinearElasticity2DExtensionPlaneStressLagrangeBasis
#ifndef NOMPIMOD
  USE MPI
#endif
  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=2
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  INTEGER(CMISSIntg),   PARAMETER ::    TIMESTEPS   = 5
  INTEGER(CMISSIntg)              ::    TIMESTEP    = 0
  REAL(CMISSRP),        PARAMETER ::    ORIGIN(3)   = [0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP]
  REAL(CMISSRP),        PARAMETER ::    LENGTH      = 120.0_CMISSRP
  REAL(CMISSRP),        PARAMETER ::    WIDTH       = 160.0_CMISSRP
  REAL(CMISSRP),        PARAMETER ::    ZERO        = 0.0_CMISSRP
  REAL(CMISSRP),        PARAMETER ::    THICKNESS   = 1.0_CMISSRP
  REAL(CMISSRP),        PARAMETER ::    EMODULE     = 30.0E6_CMISSRP
  REAL(CMISSRP),        PARAMETER ::    NU          = 0.499_CMISSRP
  REAL(CMISSRP)                   ::    BCDISP      = 0.0_CMISSRP

  !Program variables
  INTEGER(CMISSIntg)                :: NumberGlobalXElements,NumberGlobalYElements
  INTEGER(CMISSIntg)                :: MPI_IERROR
  INTEGER(CMISSIntg)                :: EquationsSetIndex
  INTEGER(CMISSIntg)                :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg)                :: node_idx,component_idx,NodeNumber,NodeDomain
  INTEGER(CMISSIntg),ALLOCATABLE    :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE    :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE    :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg)                :: LeftNormalXi,RightNormalXi,FrontNormalXi
  CHARACTER(LEN=256)                :: filename

  !CMISS variables
  TYPE(cmfe_RegionType)             :: WorldRegion
  TYPE(cmfe_CoordinateSystemType)   :: WorldCoordinateSystem
  TYPE(cmfe_BasisType)              :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType)   :: CoordinateSystem
  TYPE(cmfe_DecompositionType)      :: Decomposition
  TYPE(cmfe_EquationsType)          :: Equations
  TYPE(cmfe_EquationsSetType)       :: EquationsSet
  TYPE(cmfe_FieldType)              :: GeometricField,DependentField,MaterialField
  TYPE(cmfe_FieldsType)             :: Fields
  TYPE(cmfe_MeshType)               :: Mesh
  TYPE(cmfe_GeneratedMeshType)      :: GeneratedMesh
  TYPE(cmfe_NodesType)              :: Nodes
  TYPE(cmfe_ProblemType)            :: Problem
  TYPE(cmfe_RegionType)             :: Region
  TYPE(cmfe_FieldType)              :: EquationsSetField
  TYPE(cmfe_SolverType)             :: Solver
  TYPE(cmfe_SolverEquationsType)    :: SolverEquations

#ifdef WIN32
  ! Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  ! Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

#ifdef WIN32
  ! Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  ! Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  ! If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  ! Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."
  NumberGlobalXElements=4
  NumberGlobalYElements=4

  ! Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  ! Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  ! Coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystem,(ORIGIN),Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  ! Region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  ! Basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis, &
    & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[3,3],Err)
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  ! Generated Mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,LENGTH],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)

  ! Mesh
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  ! Decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  ! Geometry
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Undeformed",Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  ! Equations set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    &  CMFE_EQUATIONS_SET_LINEAR_ELASTICITY_TYPE, &
    & CMFE_EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  ! Dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Displacement",Err)
  DO component_idx=1,2
      CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  ! Material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,THICKNESS,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,EMODULE,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,NU,Err)

  ! Equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(EQUATIONS,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
                                              !CMFE_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                              !CMFE_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations.
  CALL cmfe_Equations_OutputTypeSet(EQUATIONS,CMFE_EQUATIONS_NO_OUTPUT,Err)
                                            !CMFE_EQUATIONS_NO_OUTPUT !<No output from the equations.
                                            !CMFE_EQUATIONS_TIMING_OUTPUT !<Timing information output.
                                            !CMFE_EQUATIONS_MATRIX_OUTPUT !<All below and equation matrices output.
                                            !CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT !<All below and Element matrices output.
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  ! Problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber, &
    & [CMFE_PROBLEM_ELASTICITY_CLASS, &
    & CMFE_PROBLEM_LINEAR_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  ! Control Loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  ! Solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(SOLVER,CMFE_SOLVER_NO_OUTPUT,Err)
                                      !CMFE_SOLVER_NO_OUTPUT !<No output from the solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                      !CMFE_SOLVER_PROGRESS_OUTPUT !<Progress output from solver routines.
                                      !CMFE_SOLVER_TIMING_OUTPUT !<Timing output from the solver routines plus below.
                                      !CMFE_SOLVER_SOLVER_OUTPUT !<Solver specific output from the solver routines plus below.
                                      !CMFE_SOLVER_MATRIX_OUTPUT !<Solver matrices output from the solver routines plus below.
!  CALL cmfe_Solver_LibraryTypeSet(SOLVER,CMFE_SOLVER_PETSC_LIBRARY,Err)
  CALL cmfe_Solver_LinearTypeSet(SOLVER,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
                                      !CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE    !<Direct linear solver type.
                                      !CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE !<Iterative linear solver type.
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  ! Solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
                                                          !CMFE_SOLVER_SPARSE_MATRICES !<Use sparse solver matrices.
                                                          !CMFE_SOLVER_FULL_MATRICES !<Use fully populated solver matrices.
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  ! Boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi, &
    & Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 0% x displacement
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    WRITE(*,*) NodeNumber
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  ! Solve multiple timesteps
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  DO TIMESTEP=1,TIMESTEPS
    WRITE(*,*) "TIME STEP ",TIMESTEP,"/",TIMESTEPS
    ! update BCs - move right surface nodes in positive x-direction
    BCDISP = 0.1*WIDTH*(TIMESTEP-1)/TIMESTEPS
    DO node_idx=1,SIZE(RightSurfaceNodes,1)
      NodeNumber=RightSurfaceNodes(node_idx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
          & 1,1,NodeNumber,1,BCDISP,Err)
      ENDIF
    ENDDO
    CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    ! The actual solve
    CALL cmfe_Problem_Solve(Problem,Err)
    ! Export solution
    IF (TIMESTEP.LE.9) THEN
      WRITE(filename, "(A30,I1)") "LinearElasticity2DPlaneStress_",TIMESTEP
      filename=trim(filename)
    ELSEIF (TIMESTEP.LE.99) THEN
      WRITE(filename, "(A30,I2)") "LinearElasticity2DPlaneStress_",TIMESTEP
      filename=trim(filename)
    ELSEIF (TIMESTEP.LE.999) THEN
      WRITE(filename, "(A30,I3)") "LinearElasticity2DPlaneStress_",TIMESTEP
      filename=trim(filename)
    ELSE
      WRITE(filename, "(A30,I4)") "LinearElasticity2DPlaneStress_",TIMESTEP
      filename=trim(filename)
    ENDIF
    CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  END DO

  ! Export final solution and topology
  CALL cmfe_Fields_ElementsExport(Fields,"LinearElasticity2DPlaneStress","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  ! Finalise
  CALL cmfe_Finalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  STOP

END PROGRAM LinearElasticity2DExtensionPlaneStressLagrangeBasis
