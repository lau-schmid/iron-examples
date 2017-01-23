
allocate(Region_parsing(1))
allocate(Region_arguments(1,num_of_Region))
Region_parsing = "REGION_ID"



allocate(CoordinateSystem_parsing(2))
allocate(CoordinateSystem_arguments(2,num_of_CoordinateSystem))
CoordinateSystem_parsing(1) = "CARTESIAN_SYSTEM_ID"
CoordinateSystem_parsing(2) = "TYPE"

allocate(ControlLoop_parsing(3))
allocate(ControlLoop_arguments(3,num_of_ControlLoop))
ControlLoop_parsing(1) = "CONTROL_LOOP_ID"
ControlLoop_parsing(2) = "TYPE"
ControlLoop_parsing(3) = "INCREMENTS"

allocate(DependentField_parsing(5))
allocate(DependentField_arguments(5,num_of_DependentField))
DependentField_parsing(1) = "DEPENDENT_FIELD_ID"
DependentField_parsing(2) = "TYPE"
DependentField_parsing(3) = "NUM_OF_COMPONENTS"
DependentField_parsing(4) = "INITIAL_VALUE_OF_STATE_VARIABLE"
DependentField_parsing(5) = "INITIAL_VALUE_OF_PRESSURE_STATE_VARIABLE"

allocate(EquationsSet_parsing(4))

EquationsSet_parsing(1) ="MODEL_ID"
EquationsSet_parsing(2) ="CLASS"
EquationsSet_parsing(3) ="TYPE"
EquationsSet_parsing(4) ="SUBTYPE"


allocate(EquationsSet_arg1(1,num_of_EquationsSet))
allocate(EquationsSet_arg2(1,num_of_EquationsSet))
allocate(EquationsSet_arg3(1,num_of_EquationsSet))
allocate(EquationsSet_arg4(1,num_of_EquationsSet))



allocate(Problem_parsing(4))
Problem_parsing(1) ="MODEL_ID"
Problem_parsing(2) ="CLASS"
Problem_parsing(3) ="TYPE"
Problem_parsing(4) ="SUBTYPE"


allocate(Problem_arg1(1,num_of_Problem))
allocate(Problem_arg2(1,num_of_Problem))
allocate(Problem_arg3(1,num_of_Problem))
allocate(Problem_arg4(1,num_of_Problem))

allocate(MaterialField_parsing(2))
MaterialField_parsing(1) ="MATERIAL_FIELD_ID"
MaterialField_parsing(2) ="MATERIAL_FIELD_PARAMETERS"

allocate(MaterialField_arg1(1,num_of_MaterialField))
allocate(MaterialField_arg2(10,num_of_MaterialField))

allocate(FiberField_parsing(2))
FiberField_parsing(1) = "FIBER_FIELD_ID"
FiberField_parsing(2) = "FIBER_FIELD_PARAMETERS"

allocate(FiberField_arg1(1,num_of_FiberField))
allocate(FiberField_arg2(10,num_of_FiberField))


allocate(Mesh_parsing(5))
Mesh_parsing(1) ="MESH_ID"
Mesh_parsing(2) ="TOPOLOGY"
Mesh_parsing(3) ="GEOMETRIC_PARAMETERS"
Mesh_parsing(4) ="NUMBER_OF_ELEMENTS"

allocate(Mesh_arg1(1,num_of_Mesh))
allocate(Mesh_arg2(1,num_of_Mesh))
allocate(Mesh_arg3(3,num_of_Mesh))
allocate(Mesh_arg4(3,num_of_Mesh))


allocate(BC_parsing(4))
BC_parsing(1) =  "BOUNDARY_CONDITIONS_ID"
BC_parsing(2) =  "DIRICHELET"
BC_parsing(3) =  "NEUMANN_CF"
BC_parsing(4) =  "NEUMANN_PRESSURE"


allocate(BC_arg1(1,num_of_BoundaryCondition))
allocate(BC_arg2(100,num_of_BoundaryCondition))
allocate(BC_arg3(100,num_of_BoundaryCondition))
allocate(BC_arg4(100,num_of_BoundaryCondition))

allocate(Solvers_arguments(8,num_of_Solver))
allocate(Solvers_parsing(8))

Solvers_parsing(1) =  "SOLVER_ID"
Solvers_parsing(2) =  "EQUATION_SET_SOLVER_TYPE"
Solvers_parsing(3) =  "PRECONDITIONER"
Solvers_parsing(4) =  "EQUATION_SOLVER_MAXITER"
Solvers_parsing(5) =  "EQUATION_SOLVER_TOLERANCE"
Solvers_parsing(6) =  "NEWTON_MAX_ITERATIONS"
Solvers_parsing(7) =  "NEWTON_TOLERANCE"
Solvers_parsing(8) =  "NONLINEAR_SOLVER"

allocate(Basis_arguments(5,num_of_Basis))
allocate(Basis_parsing(5))

Basis_parsing(1) =  "BASIS_ID"
Basis_parsing(2) =  "COMPONENTS"
Basis_parsing(3) =  "PRESSURE_BASIS "
Basis_parsing(4) =  "BASIS"
Basis_parsing(5) =  "NumberOfGaussXi"


allocate(PressureBasis_arguments(2,num_of_PressureBasis))
allocate(PressureBasis_parsing(2))

PressureBasis_parsing(1) =  "PRESSURE_BASIS_ID"
PressureBasis_parsing(2) =  "TYPE"

allocate(Output_arguments(2,1))
allocate(Output_parsing(2))

Output_parsing(1) =  "MATRIX_TYPE"
Output_parsing(2) =  "OUTPUT_TYPE"

allocate(Field_parsing(3))
Field_parsing(1) ="FIELD_ID"
Field_parsing(2) ="TYPE"
Field_parsing(3) ="VALUE"

allocate(Field_arg1(3,num_of_Field))
allocate(Field_arg2(3,num_of_Field))
allocate(Field_arg3(4,num_of_Field))


 close(12)

 open(12,file=fileplace,status="old")
 read(12,'(A)',iostat=filestat) rdline



do while (trim(rdline).NE."STOP_PARSING")

 read(12,'(A)',iostat=filestat) rdline

 call Equations_Set_parsing(EquationsSet_parsing,EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3,&
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


 call BC_parsing_subroutine(BC_parsing,BC_arg1,BC_arg2,BC_arg3,BC_arg4, rdline,num_of_dirichelet)



 call Mesh_parsing_subroutine(Mesh_parsing,Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4,rdline)



 call parsing_subroutine(PressureBasis_parsing,PressureBasis_arguments,rdline,2,"PRESSURE_BASIS")



 call parsing_subroutine(Output_parsing,Output_arguments,rdline,2,"OUTPUT")




end do


