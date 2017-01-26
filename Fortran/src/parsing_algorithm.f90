
!!!!!!!!!!!!!!!!!!!!!!!!!!! __________________ THIS HEADER FILE ALLOCATES AND INITIALIZE THE DATA STRUCTURES USED IN THE PARSING ALGORITHMS. __________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(Region_parsing(1))
allocate(Region_arguments(1,NumberOfRegion))                                                 !! Allocation of "Region BLock" related Data Structures.
Region_parsing = "REGION_ID"

allocate(CoordinateSystem_parsing(2))
allocate(CoordinateSystem_arguments(2,NumberOfCoordinateSystem))
CoordinateSystem_parsing(1) = "CARTESIAN_SYSTEM_ID"
CoordinateSystem_parsing(2) = "TYPE"

allocate(ControlLoop_parsing(3))                                                            !! Allocation of "ControlLoop BLock " related Data Structures.
allocate(ControlLoop_arguments(3,NumberOfControlLoop))
ControlLoop_parsing(1) = "CONTROL_LOOP_ID"
ControlLoop_parsing(2) = "TYPE"
ControlLoop_parsing(3) = "INCREMENTS"

allocate(DependentField_parsing(5))                                                            !! Allocation of "DependentFIeld BLock " related Data Structures.
allocate(DependentField_arguments(5,NumberOfDependentField))
DependentField_parsing(1) = "DEPENDENT_FIELD_ID"
DependentField_parsing(2) = "TYPE"
DependentField_parsing(3) = "NUM_OF_COMPONENTS"
DependentField_parsing(4) = "INITIAL_VALUE_OF_STATE_VARIABLE"
DependentField_parsing(5) = "INITIAL_VALUE_OF_PRESSURE_STATE_VARIABLE"

allocate(EquationsSet_parsing(4))                                                             !! Allocation of "EquationSet BLock " related Data Structures.
EquationsSet_parsing(1) ="MODEL_ID"
EquationsSet_parsing(2) ="CLASS"
EquationsSet_parsing(3) ="TYPE"
EquationsSet_parsing(4) ="SUBTYPE"
allocate(EquationsSet_arg1(1,NumberOfEquationsSet))
allocate(EquationsSet_arg2(1,NumberOfEquationsSet))
allocate(EquationsSet_arg3(1,NumberOfEquationsSet))
allocate(EquationsSet_arg4(1,NumberOfEquationsSet))

allocate(Problem_parsing(4))                                                                    !! Allocation of "PROBLEM"related Data Structures.
Problem_parsing(1) ="MODEL_ID"
Problem_parsing(2) ="CLASS"
Problem_parsing(3) ="TYPE"
Problem_parsing(4) ="SUBTYPE"
allocate(Problem_arg1(1,NumberOfProblem))
allocate(Problem_arg2(1,NumberOfProblem))
allocate(Problem_arg3(1,NumberOfProblem))
allocate(Problem_arg4(1,NumberOfProblem))

allocate(MaterialField_parsing(2))                                                            !! Allocation of "Material BLock " related Data Structures.
MaterialField_parsing(1) ="MATERIAL_FIELD_ID"
MaterialField_parsing(2) ="MATERIAL_FIELD_PARAMETERS"
allocate(MaterialField_arg1(1,NumberOfMaterialField))
allocate(MaterialField_arg2(10,NumberOfMaterialField))

allocate(FiberField_parsing(2))
FiberField_parsing(1) = "FIBER_FIELD_ID"                                                   !! Allocation of "Fiber FIeld BLock " related Data Structures.
FiberField_parsing(2) = "FIBER_FIELD_PARAMETERS"
allocate(FiberField_arg1(1,NumberOfFiberField))
allocate(FiberField_arg2(10,NumberOfFiberField))

allocate(Mesh_parsing(5))                                                                   !! Allocation of "Mesh BLock " related Data Structures.
Mesh_parsing(1) ="MESH_ID"
Mesh_parsing(2) ="TOPOLOGY"
Mesh_parsing(3) ="GEOMETRIC_PARAMETERS"
Mesh_parsing(4) ="NUMBER_OF_ELEMENTS"
allocate(Mesh_arg1(1,NumberOfMesh))
allocate(Mesh_arg2(1,NumberOfMesh))
allocate(Mesh_arg3(3,NumberOfMesh))
allocate(Mesh_arg4(3,NumberOfMesh))

allocate(BC_parsing(4))                                                                          !! Allocation of "BoundaryCondition BLock " related Data Structures.
BC_parsing(1) =  "BOUNDARY_CONDITIONS_ID"
BC_parsing(2) =  "DIRICHELET"
BC_parsing(3) =  "NEUMANN_CF"
BC_parsing(4) =  "NEUMANN_PRESSURE"
allocate(BC_arg1(1,NumberOfBoundaryCondition))
allocate(BC_arg2(100,NumberOfBoundaryCondition))
allocate(BC_arg3(100,NumberOfBoundaryCondition))
allocate(BC_arg4(100,NumberOfBoundaryCondition))

allocate(Solvers_arguments(8,NumberOfSolver))                                                 !! Allocation of "Solver BLock " related Data Structures.
allocate(Solvers_parsing(8))
Solvers_parsing(1) =  "SOLVER_ID"
Solvers_parsing(2) =  "EQUATION_SET_SOLVER_TYPE"
Solvers_parsing(3) =  "PRECONDITIONER"
Solvers_parsing(4) =  "EQUATION_SOLVER_MAXITER"
Solvers_parsing(5) =  "EQUATION_SOLVER_TOLERANCE"
Solvers_parsing(6) =  "NEWTON_MAX_ITERATIONS"
Solvers_parsing(7) =  "NEWTON_TOLERANCE"
Solvers_parsing(8) =  "NONLINEAR_SOLVER"

allocate(Basis_arguments(5,NumberOfBasis))                                                !! Allocation of "Basis BLock " related Data Structures.
allocate(Basis_parsing(5))

Basis_parsing(1) =  "BASIS_ID"
Basis_parsing(2) =  "COMPONENTS"
Basis_parsing(3) =  "PRESSURE_BASIS "
Basis_parsing(4) =  "BASIS"
Basis_parsing(5) =  "NumberOfGaussXi"

allocate(PressureBasis_arguments(2,NumberOfPressureBasis))                                !! Allocation of "Pressure BLock " related Data Structures.
allocate(PressureBasis_parsing(2))

PressureBasis_parsing(1) =  "PRESSURE_BASIS_ID"
PressureBasis_parsing(2) =  "TYPE"

allocate(Output_arguments(2,1))                                                                !! Allocation of "OUtput BLock " related Data Structures.
allocate(Output_parsing(2))

Output_parsing(1) =  "MATRIX_TYPE"
Output_parsing(2) =  "OUTPUT_TYPE"

allocate(Field_parsing(3))                                                                !! Allocation of "Field BLock " related Data Structures.
Field_parsing(1) ="FIELD_ID"
Field_parsing(2) ="TYPE"
Field_parsing(3) ="VALUE"

allocate(Field_arg1(3,NumberOfField))
allocate(Field_arg2(3,NumberOfField))
allocate(Field_arg3(4,NumberOfField))


 close(12)

 open(12,file=InputFile ,status="old")
 read(12,'(A)',iostat=filestat) rdline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! START PARSING THE INPUT FILE TO STORE THE INPUT ARGUMENTS FORM THE INPUT FILE IN THE RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!

do while (trim(rdline).NE."STOP_PARSING")

  read(12,'(A)',iostat=filestat) rdline

  call Equations_Set_parsing(EquationsSet_parsing,EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3, &
    & EquationsSet_arg4, rdline)

  call Problem_parsing_subroutine(Problem_parsing,Problem_arg1,Problem_arg2,Problem_arg3,Problem_arg4,rdline)

  call MaterialField_parsing_subroutine(MaterialField_parsing,MaterialField_arg1,MaterialField_arg2,rdline)

  call parsing_subroutine(DependentField_parsing,DependentField_arguments,rdline,5,"DEPENDENT_FIELD")

  call FiberField_parsing_subroutine(FiberField_parsing,FiberField_arg1,FiberField_arg2,rdline)

  call parsing_subroutine(Solvers_parsing,Solvers_arguments,rdline,8,"SOLVER_SETTINGS")

  call Field_parsing_subroutine(Field_parsing,Field_arg1,Field_arg2,Field_arg3,rdline)

  call parsing_subroutine(ControlLoop_parsing,ControlLoop_arguments,rdline,3,"CONTROL_LOOP")

  call parsing_subroutine(Basis_parsing,Basis_arguments,rdline,5,"BASIS")

  call parsing_subroutine(CoordinateSystem_parsing,CoordinateSystem_arguments,rdline,2,"COORDINATE_SYSTEM")

  call parsing_subroutine(Region_parsing,Region_arguments,rdline,1,"REGION")

  call BC_parsing_subroutine(BC_parsing,BC_arg1,BC_arg2,BC_arg3,BC_arg4, rdline,NumberOfDirichelet)

  call Mesh_parsing_subroutine(Mesh_parsing,Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4,rdline)

  call parsing_subroutine(PressureBasis_parsing,PressureBasis_arguments,rdline,2,"PRESSURE_BASIS")

  call parsing_subroutine(Output_parsing,Output_arguments,rdline,2,"OUTPUT")


end do

