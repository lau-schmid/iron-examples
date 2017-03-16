
!!!!!!!!!!!!!!!!!!!!!!!!!!! __________________ THIS HEADER FILE ALLOCATES AND INITIALIZE THE DATA STRUCTURES USED IN THE PARSING ALGORITHMS. __________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ALLOCATE(RegionKeyWords(2))
RegionKeyWords(1) = "REGION_ID"
RegionKeyWords(2) = "REGION_LABEL_SET"
!! Might need to add more functinalities later

ALLOCATE(MeshKeyWords(1))
MeshKeyWords(1) = "MESH_ID"
!! Might need to add more functinalities later

ALLOCATE(FieldsKeyWords(1))
FieldsKeyWords(1) = "FIELDS_ID"
!! Might need to add more functinalities later

ALLOCATE(GeometricFieldKeyWords(2))
GeometricFieldKeyWords(1) = "GEOMETRIC_FIELD_ID"
GeometricFieldKeyWords(2) = "GEOMETRIC_FIELD_LABEL_SET"
!! Might need to add more functinalities later


ALLOCATE(CoordinateSystemKeywords(3))
CoordinateSystemKeywords(1) = "CARTESIAN_SYSTEM_ID"
CoordinateSystemKeywords(2) = "TYPE"
CoordinateSystemKeywords(3) = "DIMENSION"

ALLOCATE(ControlLoopKeyWords(4))                                                               !! Allocation of "ControlLoop BLock " related Data Structures.
ControlLoopKeyWords(1) = "CONTROL_LOOP_ID"
ControlLoopKeyWords(2) = "TYPE"
ControlLoopKeyWords(3) = "INCREMENTS"
ControlLoopKeyWords(4) = "TIME_INCREMENTS"


ALLOCATE(DependentFieldKeyWords(6))                                                            !! Allocation of "DependentFIeld BLock " related Data Structures.
DependentFieldKeyWords(1) = "DEPENDENT_FIELD_ID"
DependentFieldKeyWords(2) = "STATE_VARIABLE"
DependentFieldKeyWords(3) = "NUM_OF_COMPONENTS"
DependentFieldKeyWords(4) = "INITIAL_VALUE_OF_STATE_VARIABLE"
DependentFieldKeyWords(5) = "INITIAL_VALUE_OF_PRESSURE_STATE_VARIABLE"
DependentFieldKeyWords(6) = "DEPENDENT_FIELD_LABEL_SET"

ALLOCATE(EquationSetKeyWords(5))                                                               !! Allocation of "EquationSet BLock " related Data Structures.
EquationSetKeyWords(1) ="MODEL_ID"
EquationSetKeyWords(2) ="CLASS"
EquationSetKeyWords(3) ="TYPE"
EquationSetKeyWords(4) ="SUBTYPE"
EquationSetKeyWords(5) ="EQUATIONS_OUTPUT_TYPE_SET"


ALLOCATE(ProblemKeywords(4))                                                                   !! Allocation of "PROBLEM"related Data Structures.
ProblemKeywords(1) ="MODEL_ID"
ProblemKeywords(2) ="CLASS"
ProblemKeywords(3) ="TYPE"
ProblemKeywords(4) ="SUBTYPE"


ALLOCATE(MaterialFieldKeyWords(3))                                                            !! Allocation of "Material BLock " related Data Structures.
MaterialFieldKeyWords(1) ="MATERIAL_FIELD_ID"
MaterialFieldKeyWords(2) ="MATERIAL_FIELD_PARAMETERS"
MaterialFieldKeyWords(3) ="MATERIAL_FIELD_LABEL_SET"

ALLOCATE(FibreFieldKeywords(3))
FibreFieldKeywords(1) = "FIBER_FIELD_ID"                                                      !! Allocation of "Fiber FIeld BLock " related Data Structures.
FibreFieldKeywords(2) = "FIBER_FIELD_PARAMETERS"
FibreFieldKeywords(3) = "FIBER_FIELD_LABEL_SET"


ALLOCATE(GeneratedMeshKeyWords(5))                                                                     !! Allocation of "Mesh BLock " related Data Structures.
GeneratedMeshKeyWords(1) ="GENERATED_MESH_ID"
GeneratedMeshKeyWords(2) ="GENERATED_MESH_TYPE"
GeneratedMeshKeyWords(3) ="GENERATED_MESH_GEOMETRIC_EXTENTS"
GeneratedMeshKeyWords(4) ="GENERATED_MESH_ORIGIN_SET"
GeneratedMeshKeyWords(5) ="GENERATED_MESH_NUMBER_OF_ELEMENTS"


ALLOCATE(BoundaryConditionKeyWords(4))                                                        !! Allocation of "BoundaryCondition BLock " related Data Structures.
BoundaryConditionKeyWords(1) =  "BOUNDARY_CONDITIONS_ID"
BoundaryConditionKeyWords(2) =  "DIRICHILET"
BoundaryConditionKeyWords(3) =  "NEUMANN_CF"
BoundaryConditionKeyWords(4) =  "NEUMANN_PRESSURE"

!! Allocation of "Solver BLock " related Data Structures.
ALLOCATE(SolverKeyWords(15))
SolverKeyWords(1) =  "SOLVER_ID"
!! solver for solving set of equations
SolverKeyWords(2) =  "LINEAR_SOLVER_TYPE"
!! direct solver
SolverKeyWords(3) = "EQUATION_DIRECT_SOLVER_TYPE"
!! properties of iterative solvers
SolverKeyWords(4) = "EQUATION_ITERATIVE_SOLVER_TYPE"
SolverKeyWords(5) = "LINEAR_ITERATIVE_RELATIVE_TOLERANCE"
SolverKeyWords(6) = "LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE"
SolverKeyWords(7) = "LINEAR_ITERATIVE_MAXIMUM_ITERATIVE"
SolverKeyWords(8) = "LINEAR_ITERATIVE_PRECONDITIONER_TYPE"
!! properties of nonlinear solver
SolverKeyWords(9) =  "NONLINEAR_SOLVER"
SolverKeyWords(10) =  "NEWTON_SOLVER_MAXIMUM_ITERATIONS"
SolverKeyWords(11) =  "NEWTON_SOLVER_ABSOLUTE_TOLERANCE"
SolverKeyWords(12) =  "NEWTON_SOLVER_RELATIVE_TOLERANCE"
!! methods for calculating jacobians
SolverKeyWords(13) =  "JACOBIAN_TYPE"
!! solver library
SolverKeyWords(14) = "SOLVER_LIBRARY_TYPE"
SolverKeyWords(15) = "SOLVER_OUTPUT_TYPE_SET"


ALLOCATE(BasisKeyWords(3))
BasisKeyWords(1) =  "BASIS_ID"
BasisKeyWords(2) =  "BASIS_INTERPOLATION_TYPE"
BasisKeyWords(3) =  "NUMBEROFGAUSSXI"


ALLOCATE(OutputKeyWords(3))                                                                  !! Allocation of "OUtput BLock " related Data Structures.
OutputKeyWords(1) =  "OUTPUT_ID"
OutputKeyWords(2) =  "NODE_EXPORT"
OutputKeywords(3) =  "ELEMENT_EXPORT"



ALLOCATE(FieldKeyWords(3))                                                                   !! Allocation of "Field BLock " related Data Structures.
FieldKeyWords(1) ="FIELD_ID"
FieldKeyWords(2) ="TYPE"
FieldKeyWords(3) ="VALUE"


ALLOCATE(DecompositionKeyWords(2))                                                           !! Allocation of "DECOMPOSITION BLock " related Data Structures.
DecompositionKeyWords(1) ="CALCULATE_ELEMENT_FACES"
DecompositionKeyWords(2) ="DECOMPOSITION_ID"

ALLOCATE(FunctionKeyWords(2))                                                           !! Allocation of "DECOMPOSITION BLock " related Data Structures.
FunctionKeyWords(1) ="FUNCTION_ID"
FunctionKeyWords(2) ="FUNCTION_CONSTANTS"

 CLOSE(12)

open(12,file=InputFile ,status="old")
READ(12,'(A)',iostat=filestat) rdline

DependentFieldId   = 1
SolverSettingsId   = 1
ControlLoopId      = 1
BasisId            = 1
CoordinateSystemId = 1
RegionId           = 1
OutputId           = 1
EquationId         = 1
ProblemId          = 1
MaterialFieldId    = 1
FibreFieldId       = 1
GeneratedMeshId     = 1
DecompositionId    = 1
SourceFieldId      = 1
GeometricFieldId   = 1
MeshId             = 1
GeneratedMeshId    = 1
FieldsId           = 1
FunctionId         = 1

ALLOCATE(BoundaryConditionId(1,NumberOfBoundaryCondition))
ALLOCATE(BoundaryDirichelet(100,NumberOfBoundaryCondition))
ALLOCATE(BoundaryTractionNeumann(100,NumberOfBoundaryCondition))
ALLOCATE(BoundaryPressureNeumann(100,NumberOfBoundaryCondition))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START PARSING THE INPUT FILE TO STORE THE INPUT ARGUMENTS FORM THE INPUT FILE IN THE RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!
DO WHILE (TRIM(rdline).NE."STOP_PARSING")

  READ(12,'(A)',iostat=filestat) rdline

  IF (NumberOfSourceField .GT. 0) THEN

   CALL GenericParsing(FieldKeyWords, "SOURCE_FIELD",3, SourceFieldId, rdline,all_SourceField(SourceFieldId)%SourceFieldId(:), &
     & all_SourceField(SourceFieldId)%SourceFieldType(:),all_SourceField(SourceFieldId)%SourceFieldComponents(:))

  END IF

  CALL GenericParsing(ControlLoopKeyWords, "CONTROL_LOOP",4, ControlLoopId, rdline, &
    & all_ControlLoop(ControlLoopId)%ControlLoopId(:), &
      & all_ControlLoop(ControlLoopId)%ControlLoopType(:),all_ControlLoop(ControlLoopId)%ControlLoopLoadIncrement(:), &
        & all_ControlLoop(ControlLoopId)%ControlLoopTimeIncrement(:))

  CALL GenericParsing(RegionKeyWords, "REGION",2, RegionId, rdline,all_Region(RegionId)%RegionId(:),  &
    all_Region(RegionId)%RegionLabel(:))

  CALL GenericParsing(FieldsKeyWords, "FIELDS",1, FieldsId, rdline,all_Fields(FieldsId)%FieldsId(:))

  CALL GenericParsing(GeometricFieldKeyWords, "GEOMETRIC_FIELD",2, GeometricFieldId, rdline, &
    & all_GeometricField(GeometricFieldId)%GeometricFieldId(:),all_GeometricField(GeometricFieldId)%GeometricFieldLabel(:))

  CALL BC_parsing_subroutine(BoundaryConditionKeyWords,BoundaryConditionID,BoundaryDirichelet,BoundaryTractionNeumann, &
    & BoundaryPressureNeumann, rdline,NumberOfDirichelet)

  CALL GenericParsing(OutputKeyWords,"OUTPUT",3,OutputId,rdline,all_Output(OutputId)%OutputId(:), &
    all_Output(OutputId)%NodeExport(:), all_Output(OutputId)%ElementExport(:) )

  CALL GenericParsing(DecompositionKeyWords,"DECOMPOSITION",2,DecompositionId,rdline,  &
    & all_Decomposition(DecompositionId)%CalculateElementFaces(:),all_Decomposition(DecompositionId)%DecompositionId(:))

  CALL GenericParsing(MeshKeyWords,"MESH",1,MeshId,rdline,all_Mesh(MeshId)%MeshId(:))

  CALL GenericParsing(EquationSetKeyWords,"EQUATIONS_SET",5,EquationId,rdline,all_EquationsSet(EquationId)%EquationSetId(:), &
    & all_EquationsSet(EquationId)%EquationSetClass(:), all_EquationsSet(EquationId)%EquationSetType(:), &
      & all_EquationsSet(EquationId)%EquationSetSubType(:),all_EquationsSet(EquationId)%EquationSetOutputTypeSet(:))

  CALL GenericParsing(ProblemKeyWords,"PROBLEM",4,ProblemId,rdline,all_Problem(ProblemId)%ProblemId(:), &
    & all_Problem(ProblemId)%ProblemClass(:), all_Problem(ProblemId)%ProblemType(:), &
      & all_Problem(ProblemId)%ProblemSubType(:))

  CALL GenericParsing(MaterialFieldKeyWords,"MATERIAL_FIELD",3,MaterialFieldId,rdline, &
    & all_MaterialField(MaterialFieldId)%MaterialFieldId(:),all_MaterialField(MaterialFieldId)%MaterialFieldParameters(:), &
      & all_MaterialField(MaterialFieldId)%MaterialFieldLabel(:))

  IF (NumberofFiberField .GT. 0) then
    CALL GenericParsing(FibreFieldKeyWords,"FIBER_FIELD",3,FibreFieldId,rdline, &
      & all_FibreField(NumberOfFiberField)%FibreFieldId(:),all_FibreField(NumberOfFiberField)%FibreFieldParameters(:), &
        & all_FibreField(FibreFieldId)%FiberFieldLabel(:))
  END IF

  CALL GenericParsing(SolverKeyWords,"SOLVER_SETTINGS",15,SolverSettingsId,rdline,all_Solver(SolverSettingsId)%SolverId(:), &
    & all_Solver(SolverSettingsId)%LinearSolverType(:),all_Solver(SolverSettingsId)%DirectSolverType(:), &
      & all_Solver(SolverSettingsId)%IterativeSolverType(:),all_Solver(SolverSettingsId)%IterativeSolverRelativeTolerance(:), &
        & all_Solver(SolverSettingsId)%IterativeSolverAbsoluteTolerance(:),  &
          & all_Solver(SolverSettingsId)%IterativeSolverMaximumIterations(:), &
            & all_Solver(SolverSettingsId)%IterativePreconditioner(:),all_Solver(SolverSettingsId)%NonlinearSolverType(:),&
              & all_Solver(SolverSettingsId)%NewtonMaximumIterations(:),all_Solver(SolverSettingsId)%NewtonAbsoluteTolerance(:),&
                & all_Solver(SolverSettingsId)%NewtonRelativeTolerance(:) , all_Solver(SolverSettingsId)%JacobianType(:), &
                  & all_Solver(SolverSettingsId)%LibraryType(:), all_Solver(SolverSettingsId)%OutputTypeSet(:))



  CALL GenericParsing(CoordinateSystemKeyWords,"COORDINATE_SYSTEM",3,CoordinateSystemId,rdline, &
    & all_CoordinateSystem(CoordinateSystemId)%CoordinateSystemId(:),  &
      & all_CoordinateSystem(CoordinateSystemId)%CoordinateSystemType(:), &
        & all_CoordinateSystem(CoordinateSystemId)%CoordinateSystemDimension(:))

  CALL GenericParsing(GeneratedMeshKeyWords,"GENERATED_MESH",5,GeneratedMeshId,rdline,  &
    & all_GeneratedMesh(GeneratedMeshId)%GeneratedMeshID(:), &
      & all_GeneratedMesh(GeneratedMeshId)%GeneratedMeshType(:), &
        & all_GeneratedMesh(GeneratedMeshId)%GeneratedMeshGeometricExtents(:), &
          & all_GeneratedMesh(GeneratedMeshId)%GeneratedMeshOriginSet(:),  &
            & all_GeneratedMesh(GeneratedMeshId)%GeneratedMeshNumberOfElements(:))

  CALL GenericParsing(DependentFieldKeyWords,"DEPENDENT_FIELD",6,DependentFieldId,rdline, &
    & all_DependentField(DependentFieldId)%DependentFieldID(:), &
      & all_DependentField(DependentFieldId)%DependentFieldStateVariable(:),  &
        & all_DependentField(DependentFieldId)%DependentFieldNumberOfComponents(:),  &
          & all_DependentField(DependentFieldId)%DependentFieldInitialValueOfStateVector(:), &
            & all_DependentField(DependentFieldId)%DependentFieldInitialValueOfStateScalar(:), &
              & all_DependentField(DependentFieldId)%DependentFieldLabel(:))


  CALL GenericParsing(BasisKeyWords,"BASIS",3,BasisId,rdline,all_Basis(BasisId)%BasisId(:),  &
    &  all_Basis(BasisId)%BasisInterpolationType(:), all_Basis(BasisId)%BasisNumberOfGaussPoints(:))

   CALL GenericParsing(FunctionKeyWords,"FUNCTION",2,FunctionId,rdline,all_Function(FunctionId)%FunctionId(:),  &
    & all_Function(FunctionId)%FunctionConstants(:) )

END DO

