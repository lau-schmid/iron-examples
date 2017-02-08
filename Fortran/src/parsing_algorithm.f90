
!!!!!!!!!!!!!!!!!!!!!!!!!!! __________________ THIS HEADER FILE ALLOCATES AND INITIALIZE THE DATA STRUCTURES USED IN THE PARSING ALGORITHMS. __________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ALLOCATE(RegionKeyWords(1))
RegionKeyWords(1) = "REGION_ID"
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


ALLOCATE(DependentFieldKeyWords(5))                                                            !! Allocation of "DependentFIeld BLock " related Data Structures.
DependentFieldKeyWords(1) = "Dependent_FIELD_ID"
DependentFieldKeyWords(2) = "STATE_VARIABLE"
DependentFieldKeyWords(3) = "NUM_OF_COMPONENTS"
DependentFieldKeyWords(4) = "INITIAL_VALUE_OF_STATE_VARIABLE"
DependentFieldKeyWords(5) = "INITIAL_VALUE_OF_PRESSURE_STATE_VARIABLE"


ALLOCATE(EquationSetKeyWords(4))                                                               !! Allocation of "EquationSet BLock " related Data Structures.
EquationSetKeyWords(1) ="MODEL_ID"
EquationSetKeyWords(2) ="CLASS"
EquationSetKeyWords(3) ="TYPE"
EquationSetKeyWords(4) ="SUBTYPE"


ALLOCATE(ProblemKeywords(4))                                                                   !! Allocation of "PROBLEM"related Data Structures.
ProblemKeywords(1) ="MODEL_ID"
ProblemKeywords(2) ="CLASS"
ProblemKeywords(3) ="TYPE"
ProblemKeywords(4) ="SUBTYPE"


ALLOCATE(MaterialFieldKeyWords(2))                                                            !! Allocation of "Material BLock " related Data Structures.
MaterialFieldKeyWords(1) ="MATERIAL_FIELD_ID"
MaterialFieldKeyWords(2) ="MATERIAL_FIELD_PARAMETERS"


ALLOCATE(FibreFieldKeywords(2))
FibreFieldKeywords(1) = "FIBER_FIELD_ID"                                                      !! Allocation of "Fiber FIeld BLock " related Data Structures.
FibreFieldKeywords(2) = "FIBER_FIELD_PARAMETERS"


ALLOCATE(MeshKeyWords(5))                                                                     !! Allocation of "Mesh BLock " related Data Structures.
MeshKeyWords(1) ="MESH_ID"
MeshKeyWords(2) ="TOPOLOGY"
MeshKeyWords(3) ="GEOMETRIC_PARAMETERS"
MeshKeyWords(4) ="NUMBER_OF_ELEMENTS"


ALLOCATE(BoundaryConditionKeyWords(4))                                                        !! Allocation of "BoundaryCondition BLock " related Data Structures.
BoundaryConditionKeyWords(1) =  "BOUNDARY_CONDITIONS_ID"
BoundaryConditionKeyWords(2) =  "DIRICHELET"
BoundaryConditionKeyWords(3) =  "NEUMANN_CF"
BoundaryConditionKeyWords(4) =  "NEUMANN_PRESSURE"


ALLOCATE(SolverKeyWords(9))                                                                   !! Allocation of "Solver BLock " related Data Structures.
SolverKeyWords(1) =  "SOLVER_ID"
SolverKeyWords(2) =  "EQUATION_SET_SOLVER_TYPE"
SolverKeyWords(3) =  "PRECONDITIONER"
SolverKeyWords(4) =  "EQUATION_SOLVER_MAXITER"
SolverKeyWords(5) =  "EQUATION_SOLVER_TOLERANCE"
SolverKeyWords(6) =  "NEWTON_MAX_ITERATIONS"
SolverKeyWords(7) =  "NEWTON_TOLERANCE"
SolverKeyWords(8) =  "NONLINEAR_SOLVER"
SolverKeyWords(9) =  "JACOBIAN_TYPE"


ALLOCATE(BasisKeyWords(4))
BasisKeyWords(1) =  "BASIS_ID"
BasisKeyWords(2) =  "COMPONENTS"
BasisKeyWords(3) =  "BASIS"
BasisKeyWords(4) =  "NumberOfGaussXi"


ALLOCATE(OutputKeyWords(2))                                                                  !! Allocation of "OUtput BLock " related Data Structures.
OutputKeyWords(1) =  "MATRIX_TYPE"
OutputKeyWords(2) =  "OUTPUT_TYPE"

ALLOCATE(FieldKeyWords(3))                                                                   !! Allocation of "Field BLock " related Data Structures.
FieldKeyWords(1) ="FIELD_ID"
FieldKeyWords(2) ="TYPE"
FieldKeyWords(3) ="VALUE"


ALLOCATE(DecompositionKeyWords(3))                                                           !! Allocation of "DECOMPOSITION BLock " related Data Structures.
DecompositionKeyWords(1) ="NUMBER_OF_DOMAINS"
DecompositionKeyWords(2) ="FACE_ACTIVE"

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
CoordinateSystemId = 1
MeshSettingsId     = 1
DecompositionId    = 1
SourceFieldId      = 1

ALLOCATE(Output%MatrixType(1,NumberOfOutput))
ALLOCATE(Output%OutputType(1,NumberOfOutput))

ALLOCATE(EquationsSet%EquationSetId(1,NumberOfEquationsSet))
ALLOCATE(EquationsSet%EquationSetClass(1,NumberOfEquationsSet))
ALLOCATE(EquationsSet%EquationSetType(1,NumberOfEquationsSet))
ALLOCATE(EquationsSet%EquationSetSubType(1,NumberOfEquationsSet))

ALLOCATE(Problem%ProblemId(1,NumberOfProblem))
ALLOCATE(Problem%ProblemClass(1,NumberOfProblem))
ALLOCATE(Problem%ProblemType(1,NumberOfProblem))
ALLOCATE(Problem%ProblemSubType(1,NumberOfProblem))

ALLOCATE(MaterialField%MaterialFieldId(1,NumberOfMaterialField))
ALLOCATE(MaterialField%MaterialFieldParameters(10,NumberOfMaterialField))

ALLOCATE(FibreField%FibreFieldId(1,NumberOfFiberField))
ALLOCATE(FibreField%FibreFieldParameters(3,NumberOfFiberField))

ALLOCATE(Solver%SolverId(1,NumberOfSolver))
ALLOCATE(Solver%Preconditioner(1,NumberOfSolver))
ALLOCATE(Solver%SolverType(1,NumberOfSolver))
ALLOCATE(Solver%EquationSolverMaximumIteration(1,NumberOfSolver))

ALLOCATE(Solver%EquationSolverTolerance(1,NumberOfSolver))
ALLOCATE(Solver%NewtonMaxIterations(1,NumberOfSolver))
ALLOCATE(Solver%NewtonTolerance(1,NumberOfSolver))
ALLOCATE(Solver%NewtonJacobianType(1,NumberOfSolver))
ALLOCATE(Solver%NonlinearSolver(1,NumberOfSolver))

ALLOCATE(CoordinateSystem%CoordinateSystemId(1,NumberOfCoordinateSystem))
ALLOCATE(CoordinateSystem%CoordinateSystemType(1,NumberOfCoordinateSystem))
ALLOCATE(CoordinateSystem%CoordinateSystemDimension(1,NumberOfCoordinateSystem))

ALLOCATE(Mesh%MeshID(1,NumberOfMesh))
ALLOCATE(Mesh%MeshTopology(1,NumberOfMesh))
ALLOCATE(Mesh%MeshGeometricParameters(3,NumberOfMesh))
ALLOCATE(Mesh%MeshNumberOfElements(3,NumberOfMesh))

ALLOCATE(Decomposition%NumberOfDomains(1,1))         !! Always 1
ALLOCATE(Decomposition%DecompositionFaceActive(1,1)) !! Always 1

ALLOCATE(DependentField%DependentFieldID(1,NumberOfDependentField))
ALLOCATE(DependentField%DependentFieldStateVariable(1,NumberOfDependentField))
ALLOCATE(DependentField%DependentFieldNumberOfComponents(1,NumberOfDependentField))
ALLOCATE(DependentField%DependentFieldInitialValueOfStateVector(3,NumberOfDependentField))
ALLOCATE(DependentField%DependentFieldInitialValueOfStateScalar(1,NumberOfDependentField))


ALLOCATE(Basis%BasisId(1,NumberOfBasis))
ALLOCATE(Basis%BasisType(1,NumberOfBasis))
ALLOCATE(Basis%BasisComponents(1,NumberOfBasis))
ALLOCATE(Basis%BasisNumberOfGaussPoints(3,NumberOfBasis))

ALLOCATE(Region%RegionId(1,NumberOfRegion))

ALLOCATE(ControlLoop%ControlLoopId(1,NumberOfControlLoop))
ALLOCATE(ControlLoop%ControlLoopType(1,NumberOfControlLoop))
ALLOCATE(ControlLoop%ControlLoopTimeIncrement(3,NumberOfControlLoop))
ALLOCATE(ControlLoop%ControlLoopLoadIncrement(3,NumberOfControlLoop))

ALLOCATE(SourceField%SourceFieldId(1,NumberOfSourceField))
ALLOCATE(SourceField%SourceFieldType(1,NumberOfSourceField))
ALLOCATE(SourceField%SourceFieldComponents(3,NumberOfSourceField))

ALLOCATE(BoundaryConditionId(1,NumberOfBoundaryCondition))
ALLOCATE(BoundaryDirichelet(100,NumberOfBoundaryCondition))
ALLOCATE(BoundaryTractionNeumann(100,NumberOfBoundaryCondition))
ALLOCATE(BoundaryPressureNeumann(100,NumberOfBoundaryCondition))

! PLease remove it later

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START PARSING THE INPUT FILE TO STORE THE INPUT ARGUMENTS FORM THE INPUT FILE IN THE RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!
!Solver%SolverType(:),Solver%Preconditioner(:),Solver%EquationSolverTolerance(:),Solver%NewtonMaxIterations(:), &
!                                                  Solver%NewtonTolerance(:),Solver%JacobianType(:)
DO WHILE (TRIM(rdline).NE."STOP_PARSING")

  READ(12,'(A)',iostat=filestat) rdline

  CALL GenericParsing(FieldKeyWords, "SOURCE_FIELD",3, SourceFieldId, rdline,SourceField%SourceFieldId(:,SourceFieldId), &
                        & SourceField%SourceFieldType(:,SourceFieldId),SourceField%SourceFieldComponents(:,SourceFieldId))

  CALL GenericParsing(ControlLoopKeyWords, "CONTROL_LOOP",4, ControlLoopId, rdline,ControlLoop%ControlLoopId(:,ControlLoopId), &
                        & ControlLoop%ControlLoopType(:,ControlLoopId),ControlLoop%ControlLoopLoadIncrement(:,ControlLoopId), &
                          & ControlLoop%ControlLoopTimeIncrement(:,ControlLoopId))

  CALL GenericParsing(RegionKeyWords, "REGION",1, RegionId, rdline,Region%RegionId(:,RegionId))

  CALL BC_parsing_subroutine(BoundaryConditionKeyWords,BoundaryConditionID,BoundaryDirichelet,BoundaryTractionNeumann, &
                               & BoundaryPressureNeumann, rdline,NumberOfDirichelet)

  CALL GenericParsing(OutputKeyWords,"OUTPUT",2,OutputId,rdline,Output%MatrixType(:,OutputId), &
                        & Output%OutputType(:,OutputId))

  CALL GenericParsing(EquationSetKeyWords,"DECOMPOSITION",2,DecompositionId,rdline, &
                       & Decomposition%NumberOfDomains(:,DecompositionId), &
                        & Decomposition%DecompositionFaceActive(:,DecompositionId))

  CALL GenericParsing(EquationSetKeyWords,"EQUATIONS_SET",4,EquationId,rdline,EquationsSet%EquationSetId(:,EquationId), &
                       & EquationsSet%EquationSetClass(:,EquationId), EquationsSet%EquationSetType(:,EquationId), &
                         & EquationsSet%EquationSetSubType(:,EquationId))

  CALL GenericParsing(ProblemKeyWords,"PROBLEM",4,ProblemId,rdline,Problem%ProblemId(:,ProblemId), &
                       & Problem%ProblemClass(:,ProblemId), Problem%ProblemType(:,ProblemId), Problem%ProblemSubType(:,ProblemId))

  CALL GenericParsing(MaterialFieldKeyWords,"MATERIAL_FIELD",2,MaterialFieldId,rdline, &
                       & MaterialField%MaterialFieldId(:,MaterialFieldId),MaterialField%MaterialFieldParameters(:,MaterialFieldId))

  CALL GenericParsing(FibreFieldKeyWords,"FIBER_FIELD",2,FibreFieldId,rdline, &
                       & FibreField%FibreFieldId(:,FibreFieldId),FibreField%FibreFieldParameters(:,FibreFieldId))

  CALL GenericParsing(SolverKeyWords,"SOLVER_SETTINGS",9,SolverSettingsId,rdline,Solver%SolverId(:,SolverSettingsId), &
                        & Solver%SolverType(:,SolverSettingsId), &
                          & Solver%Preconditioner(:,SolverSettingsId), &
                            & Solver%EquationSolverMaximumIteration(:,SolverSettingsId), &
                              & Solver%EquationSolverTolerance(:,SolverSettingsId), &
                                & Solver%NewtonMaxIterations(:,SolverSettingsId),  &
                                  & Solver%NewtonTolerance(:,SolverSettingsId), &
                                    & Solver%NonlinearSolver(:,SolverSettingsId),Solver%NewtonJacobianType(:,SolverSettingsId))

  CALL GenericParsing(CoordinateSystemKeyWords,"COORDINATE_SYSTEM",3,CoordinateSystemId,rdline, &
                       & CoordinateSystem%CoordinateSystemId(:,CoordinateSystemId), &
                         & CoordinateSystem%CoordinateSystemType(:,CoordinateSystemId), &
                           & CoordinateSystem%CoordinateSystemDimension(:,CoordinateSystemId))

  CALL GenericParsing(MeshKeyWords,"MESH",4,MeshSettingsId,rdline, Mesh%MeshId(:,MeshSettingsId), &
                        & Mesh%MeshTopology(:,MeshSettingsId), Mesh%MeshGeometricParameters(:,MeshSettingsId), &
                          & Mesh%MeshNumberOfElements(:,MeshSettingsId))

  CALL GenericParsing(DependentFieldKeyWords,"DEPENDENT_FIELD",5,DependentFieldId,rdline, &
                        & DependentField%DependentFieldID(:,DependentFieldId), &
                          & DependentField%DependentFieldStateVariable(:,DependentFieldId),  &
                            & DependentField%DependentFieldNumberOfComponents(:,DependentFieldId),  &
                              & DependentField%DependentFieldInitialValueOfStateVector(:,DependentFieldId), &
                                & DependentField%DependentFieldInitialValueOfStateScalar(:,DependentFieldId))


  CALL GenericParsing(BasisKeyWords,"BASIS",4,BasisId,rdline,Basis%BasisId(:,BasisId), Basis%BasisComponents(:,BasisId), &
                       & Basis%BasisType(BasisId,:),Basis%BasisNumberOfGaussPoints(:,BasisId))

END DO

