
MODULE STIMULATION_BOUNDARY_CONDITIONS

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
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

CONTAINS


SUBROUTINE SetBoundaryConditions()
  INTEGER(CMISSIntg) :: NodeMUserNumber
  INTEGER(CMISSIntg) :: NodeIdx

  !Prescribe boundary conditions for monodomain
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsM,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsM,BoundaryConditionsM,Err)
  
  !PRINT*, "Print SolverEquationsM Solver"
  
  !CALL cmfe_PrintSolver(SolverEquationsM, 5, 10, Err)
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  
  !SolverEquationsM
  !SolverEquationsM%solverEquations%SOLVER%SOLVERS%SOLVERS(2)% &
  !  & PTR%LINKED_SOLVERS(1)%PTR%LINEAR_SOLVER%ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE = 1e-5 


  !PRINT*, "After changing tolerance"
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  ! Finalize solver, set PETSC parameters, such as tolerances, solver type etc.
  ! The parameters are already set in the data structure SolverEquations
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsM,Err)

  !Prescribe boundary conditions for Finite Elasticity (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsFE,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsFE,BoundaryConditionsFE,Err)

  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMeshFE,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)

  !Set x=0 nodes to no x displacment in x.
  DO NodeIdx = 1, SIZE(LeftSurfaceNodes, 1)
    NodeNumber = LeftSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE, NodeNumber, 1, NodeDomain, Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      !                                    BOUNDARY_CONDITIONS,  FIELD,            VARIABLE_TYPE,
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE, DependentFieldFE, CMFE_FIELD_U_VARIABLE_TYPE, &
      !   VERSION_NUMBER, DERIVATIVE_NUMBER, USER_NODE_NUMBER, COMPONENT_NUMBER, CONDITION,
        & 1,              1,                 NodeNumber,       1,                CMFE_BOUNDARY_CONDITION_FIXED, &   ! CMFE_BOUNDARY_CONDITION_FIXED_USER_CONTROLLED
      ! VALUE
        & 0.0_CMISSRP, Err)
    ENDIF
  ENDDO

 !Set x=PhysicalWidth nodes to InitialStretch x displacement
  DO NodeIdx = 1, SIZE(RightSurfaceNodes, 1)
    NodeNumber = RightSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE, NodeNumber, 1, NodeDomain, Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED, PhysicalLength*InitialStretch, Err)
        
      !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber,&
      !  & 1,CMFE_BOUNDARY_CONDITION_FIXED, 0.0_CMISSRP, Err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO NodeIdx = 1, SIZE(FrontSurfaceNodes, 1)
    NodeNumber = FrontSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO NodeIdx = 1, SIZE(BottomSurfaceNodes, 1)
    NodeNumber = BottomSurfaceNodes(NodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(DecompositionFE,NodeNumber,1,NodeDomain,Err)
    IF (NodeDomain == ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsFE,DependentFieldFE,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsFE,Err)

END SUBROUTINE SetBoundaryConditions


SUBROUTINE SetStimulationAtNodes(StimValuePerNode)
  REAL(CMISSRP), INTENT(IN) :: StimValuePerNode
  INTEGER(CMISSIntg) :: StimulatedNodeBegin, StimulatedNodeEnd, StimulatedNodeNo
  LOGICAL :: CurrentFibreFires
  
  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    !>Find the component ID in the given field for the variable defined by the given variable ID in the provided CellML environment.
    !! This generic routine will be used to map variable ID's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
    !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
    !! - is the model ID also needed?
    !! - because the CellML fields should all be set up to allow direct use in the CellML code, the component number matches the index of the given variable in its associated array in the CellML generated code.
    !SUBROUTINE CELLML_FIELD_COMPONENT_GET_C(CELLML,MODEL_INDEX,CELLML_FIELD_TYPE,VARIABLE_ID,COMPONENT_USER_NUMBER,ERR,ERROR,*)
      !Argument variables
      !TYPE(CELLML_TYPE), POINTER :: CELLML !<The CellML environment object from which to get the field component.
      !INTEGER(INTG), INTENT(IN) :: MODEL_INDEX !<The index of the CellML model to map from.
      !INTEGER(INTG), INTENT(IN) :: CELLML_FIELD_TYPE !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CMISS_CELLML
      !CHARACTER(LEN=*), INTENT(IN) :: VARIABLE_ID !<The ID of the model variable which needs to be located in the provided field.
      !INTEGER(INTG), INTENT(OUT) :: COMPONENT_USER_NUMBER !<On return, the field component for the model variable defined by the given ID.
    CALL cmfe_CellML_FieldComponentGet(CellMLEnvironment,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
      & "wal_environment/I_HH",StimComponent,Err)
  
  ELSEIF (ModelType == 1) THEN ! 3, "MultiPhysStrain", numerically more stable
    CALL cmfe_CellML_FieldComponentGet(CellMLEnvironment,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
      & "Aliev_Panfilov/I_HH",StimComponent,Err)
      
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    CALL cmfe_CellML_FieldComponentGet(CellMLEnvironment,shortenModelIndex,CMFE_CELLML_PARAMETERS_FIELD, &
      & "Aliev_Panfilov/I_HH",StimComponent,Err)
  ENDIF
  
  NumberFiringFibres = 0
  !loop over all neuromuscular junctions (middle point of the fibres)
  DO FibreNo = 1, NumberOfFibres
    CurrentFibreFires = .FALSE.
    
    ! get middle point of fibre
    JunctionNodeNo = (FibreNo-1) * NumberOfNodesPerLongFibre + (NumberOfNodesPerLongFibre+1)/2
    
    ! add innervation zone offset
    JunctionNodeNo = JunctionNodeNo + InnervationZoneOffset(FibreNo)
    
    ! compute first node to stimulation
    StimulatedNodeBegin = JunctionNodeNo - NumberStimulatedNodesPerFibre/2
    
    ! assert limits
    StimulatedNodeBegin = MIN(FibreNo*NumberOfNodesPerLongFibre-NumberStimulatedNodesPerFibre, &
      & MAX((FibreNo-1)*NumberOfNodesPerLongFibre, StimulatedNodeBegin))
    
    StimulatedNodeEnd = StimulatedNodeBegin + NumberStimulatedNodesPerFibre-1
    
    ! loop over nodes of current fibre to be stimulated
    DO StimulatedNodeNo = StimulatedNodeBegin, StimulatedNodeEnd
      
      !                                     decomposition,  nodeUserNumber, meshComponentNumber, domain
      CALL cmfe_Decomposition_NodeDomainGet(DecompositionM, StimulatedNodeNo, 1,                   NodeDomain, Err)
      IF (NodeDomain == ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetGetNode(IndependentFieldM,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & StimulatedNodeNo,1,MotorUnitRank,Err)

        IF ((MotorUnitRank <= 0) .OR. (MotorUnitRank >= 101)) THEN
          PRINT*, ComputationalNodeNumber,": Warning! MotorUnitRank=",MotorUnitRank,", set to 100. Fibre ", FibreNo, &
             & ", Node", StimulatedNodeNo
          MotorUnitRank=100
        ELSE
          !PRINT*, "MotorUnitFiringTimes row k=", k, ": MU rank=", MotorUnitRank, ", StimComponent=",StimComponent
          MotorUnitFires = MotorUnitFiringTimes(k, MotorUnitRank)   ! determine if mu fires
          IF (MotorUnitFires == 1) THEN
            !PRINT*, "Fibre ",FibreNo,": MU ", MotorUnitRank, " fires, StimulatedNodeNo=",StimulatedNodeNo, &
            !  & ", StimComponent=",StimComponent,", StimValue=", StimValue

            CurrentFibreFires = .TRUE.
            CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField, &
              & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,StimulatedNodeNo,StimComponent,StimValuePerNode,Err)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
    ! If there was a stimulated node in the current fibre, increase number of stimulated fibres counter for output
    IF (CurrentFibreFires) NumberFiringFibres = NumberFiringFibres + 1
  ENDDO

END SUBROUTINE SetStimulationAtNodes
  

END MODULE STIMULATION_BOUNDARY_CONDITIONS