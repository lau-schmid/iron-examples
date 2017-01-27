
!!!!!!!!!!!!!!!!!!!!!!!!!!! __________________ THIS HEADER FILE ALLOCATES AND INITIALIZE THE DATA STRUCTURES USED IN THE PARSING ALGORITHMS. __________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ALLOCATE(Region_parsing(1))
ALLOCATE(Region_arguments(1,NumberOfRegion))                                                 !! Allocation of "Region BLock" related Data Structures.
Region_parsing = "REGION_ID"

ALLOCATE(CoordinateSystem_parsing(2))
ALLOCATE(CoordinateSystem_arguments(2,NumberOfCoordinateSystem))
CoordinateSystem_parsing(1) = "CARTESIAN_SYSTEM_ID"
CoordinateSystem_parsing(2) = "TYPE"

ALLOCATE(ControlLoop_parsing(3))                                                            !! Allocation of "ControlLoop BLock " related Data Structures.
ALLOCATE(ControlLoop_arguments(3,NumberOfControlLoop))
ControlLoop_parsing(1) = "CONTROL_LOOP_ID"
ControlLoop_parsing(2) = "TYPE"
ControlLoop_parsing(3) = "INCREMENTS"

ALLOCATE(DependentField_parsing(5))                                                            !! Allocation of "DependentFIeld BLock " related Data Structures.
ALLOCATE(DependentField_arguments(5,NumberOfDependentField))
DependentField_parsing(1) = "Dependent_FIELD_ID"
DependentField_parsing(2) = "STATE_VARIABLE"
DependentField_parsing(3) = "NUM_OF_COMPONENTS"
DependentField_parsing(4) = "INITIAL_VALUE_OF_STATE_VARIABLE"
DependentField_parsing(5) = "INITIAL_VALUE_OF_PRESSURE_STATE_VARIABLE"


ALLOCATE(EquationsSet_parsing(4))                                                             !! Allocation of "EquationSet BLock " related Data Structures.
EquationsSet_parsing(1) ="MODEL_ID"
EquationsSet_parsing(2) ="CLASS"
EquationsSet_parsing(3) ="TYPE"
EquationsSet_parsing(4) ="SUBTYPE"
ALLOCATE(EquationsSet_arg1(1,NumberOfEquationsSet))
ALLOCATE(EquationsSet_arg2(1,NumberOfEquationsSet))
ALLOCATE(EquationsSet_arg3(1,NumberOfEquationsSet))
ALLOCATE(EquationsSet_arg4(1,NumberOfEquationsSet))

ALLOCATE(Problem_parsing(4))                                                                    !! Allocation of "PROBLEM"related Data Structures.
Problem_parsing(1) ="MODEL_ID"
Problem_parsing(2) ="CLASS"
Problem_parsing(3) ="TYPE"
Problem_parsing(4) ="SUBTYPE"
ALLOCATE(Problem_arg1(1,NumberOfProblem))
ALLOCATE(Problem_arg2(1,NumberOfProblem))
ALLOCATE(Problem_arg3(1,NumberOfProblem))
ALLOCATE(Problem_arg4(1,NumberOfProblem))

ALLOCATE(MaterialField_parsing(2))                                                            !! Allocation of "Material BLock " related Data Structures.
MaterialField_parsing(1) ="MATERIAL_FIELD_ID"
MaterialField_parsing(2) ="MATERIAL_FIELD_PARAMETERS"
ALLOCATE(MaterialField_arg1(1,NumberOfMaterialField))
ALLOCATE(MaterialField_arg2(10,NumberOfMaterialField))

ALLOCATE(FiberField_parsing(2))
FiberField_parsing(1) = "FIBER_FIELD_ID"                                                   !! Allocation of "Fiber FIeld BLock " related Data Structures.
FiberField_parsing(2) = "FIBER_FIELD_PARAMETERS"
ALLOCATE(FiberField_arg1(1,NumberOfFiberField))
ALLOCATE(FiberField_arg2(10,NumberOfFiberField))

ALLOCATE(Mesh_parsing(5))                                                                   !! Allocation of "Mesh BLock " related Data Structures.
Mesh_parsing(1) ="MESH_ID"
Mesh_parsing(2) ="TOPOLOGY"
Mesh_parsing(3) ="GEOMETRIC_PARAMETERS"
Mesh_parsing(4) ="NUMBER_OF_ELEMENTS"
ALLOCATE(Mesh_arg1(1,NumberOfMesh))
ALLOCATE(Mesh_arg2(1,NumberOfMesh))
ALLOCATE(Mesh_arg3(3,NumberOfMesh))
ALLOCATE(Mesh_arg4(3,NumberOfMesh))

ALLOCATE(BC_parsing(4))                                                                          !! Allocation of "BoundaryCondition BLock " related Data Structures.
BC_parsing(1) =  "BOUNDARY_CONDITIONS_ID"
BC_parsing(2) =  "DIRICHELET"
BC_parsing(3) =  "NEUMANN_CF"
BC_parsing(4) =  "NEUMANN_PRESSURE"
ALLOCATE(BC_arg1(1,NumberOfBoundaryCondition))
ALLOCATE(BC_arg2(100,NumberOfBoundaryCondition))
ALLOCATE(BC_arg3(100,NumberOfBoundaryCondition))
ALLOCATE(BC_arg4(100,NumberOfBoundaryCondition))

ALLOCATE(Solvers_arguments(9,NumberOfSolver))                                                 !! Allocation of "Solver BLock " related Data Structures.
ALLOCATE(Solvers_parsing(9))
Solvers_parsing(1) =  "SOLVER_ID"
Solvers_parsing(2) =  "EQUATION_SET_SOLVER_TYPE"
Solvers_parsing(3) =  "PRECONDITIONER"
Solvers_parsing(4) =  "EQUATION_SOLVER_MAXITER"
Solvers_parsing(5) =  "EQUATION_SOLVER_TOLERANCE"
Solvers_parsing(6) =  "NEWTON_MAX_ITERATIONS"
Solvers_parsing(7) =  "NEWTON_TOLERANCE"
Solvers_parsing(8) =  "NONLINEAR_SOLVER"
Solvers_parsing(9) =  "JACOBIAN_TYPE"

ALLOCATE(Basis_arguments(5,NumberOfBasis))                                                !! Allocation of "Basis BLock " related Data Structures.
ALLOCATE(Basis_parsing(5))

Basis_parsing(1) =  "BASIS_ID"
Basis_parsing(2) =  "COMPONENTS"
Basis_parsing(3) =  "PRESSURE_BASIS"
Basis_parsing(4) =  "BASIS"
Basis_parsing(5) =  "NumberOfGaussXi"

ALLOCATE(PressureBasis_arguments(2,NumberOfPressureBasis))                                !! Allocation of "Pressure BLock " related Data Structures.
ALLOCATE(PressureBasis_parsing(2))

PressureBasis_parsing(1) =  "PRESSURE_BASIS_ID"
PressureBasis_parsing(2) =  "TYPE"

ALLOCATE(Output_arguments(2,1))                                                                !! Allocation of "OUtput BLock " related Data Structures.
ALLOCATE(Output_parsing(2))

Output_parsing(1) =  "MATRIX_TYPE"
Output_parsing(2) =  "OUTPUT_TYPE"

ALLOCATE(Field_parsing(3))                                                                !! Allocation of "Field BLock " related Data Structures.
Field_parsing(1) ="FIELD_ID"
Field_parsing(2) ="TYPE"
Field_parsing(3) ="VALUE"

ALLOCATE(Field_arg1(3,NumberOfField))
ALLOCATE(Field_arg2(3,NumberOfField))
ALLOCATE(Field_arg3(4,NumberOfField))

ALLOCATE(DecompositionParsing(3))                                                                !! Allocation of "DECOMPOSITION BLock " related Data Structures.
DecompositionParsing(1) ="NUMBER_OF_DOMAINS"
DecompositionParsing(2) ="TYPE"

ALLOCATE(DecompositionArguments(3,1))



 CLOSE(12)

open(12,file=InputFile ,status="old")
READ(12,'(A)',iostat=filestat) rdline

DependentFieldId = 0
SolverId = 0
ControlLoopId = 0
BasisId = 0
CoordinateSystemId = 0
RegionId = 0
OutputId = 0





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START PARSING THE INPUT FILE TO STORE THE INPUT ARGUMENTS FORM THE INPUT FILE IN THE RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!

DO WHILE (TRIM(rdline).NE."STOP_PARSING")

  READ(12,'(A)',iostat=filestat) rdline

  CALL Equations_Set_parsing(EquationsSet_parsing,EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3, &
    & EquationsSet_arg4, rdline)

  CALL Problem_parsing_subroutine(Problem_parsing,Problem_arg1,Problem_arg2,Problem_arg3,Problem_arg4,rdline)

  CALL MaterialField_parsing_subroutine(MaterialField_parsing,MaterialField_arg1,MaterialField_arg2,rdline)

  CALL parsing_subroutine(DependentField_parsing,DependentField_arguments,rdline,5,"DEPENDENT_FIELD",DependentFieldId)

  CALL FiberField_parsing_subroutine(FiberField_parsing,FiberField_arg1,FiberField_arg2,rdline)

  CALL parsing_subroutine(Solvers_parsing,Solvers_arguments,rdline,9,"SOLVER_SETTINGS",SolverId)

  CALL Field_parsing_subroutine(Field_parsing,Field_arg1,Field_arg2,Field_arg3,rdline)

  CALL parsing_subroutine(ControlLoop_parsing,ControlLoop_arguments,rdline,3,"CONTROL_LOOP",ControlLoopId)

  CALL parsing_subroutine(Basis_parsing,Basis_arguments,rdline,6,"BASIS",BasisId)

  CALL parsing_subroutine(CoordinateSystem_parsing,CoordinateSystem_arguments,rdline,2,"COORDINATE_SYSTEM",CoordinateSystemId)

  CALL parsing_subroutine(Region_parsing,Region_arguments,rdline,1,"REGION",RegionId)

  CALL BC_parsing_subroutine(BC_parsing,BC_arg1,BC_arg2,BC_arg3,BC_arg4, rdline,NumberOfDirichelet)

  CALL Mesh_parsing_subroutine(Mesh_parsing,Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4,rdline)

  CALL parsing_subroutine(PressureBasis_parsing,PressureBasis_arguments,rdline,6,"PRESSURE_BASIS",BasisId)

  CALL parsing_subroutine(Output_parsing,Output_arguments,rdline,2,"OUTPUT",OutputId)

  CALL parsing_subroutine(Output_parsing,Output_arguments,rdline,2,"DECOMPOSTION",DecompositionId)

END DO

