
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! EQUATION_SET block defines the governing equation!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! It creates an object of derived type CMFE_EquationsSetType !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_EQUATIONS_SET

MODEL_ID
FLUID

EQUATIONS_SET_CLASS
EQUATIONS_SET_FLUID_MECHANICS_CLASS

EQUATIONS_SET_TYPE
EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE


EQUATIONS_SET_SUBTYPE
EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE


EQUATIONS_OUTPUT_TYPE_SET
EQUATIONS_NO_OUTPUT


END_EQUATIONS_SET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! PROBLEM block defines the governing equation!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! It creates an object of derived type CMFE_ProblemType !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_PROBLEM

PROBLEM_ID
FLUID

PROBLEM_CLASS
PROBLEM_FLUID_MECHANICS_CLASS

PROBLEM_TYPE
PROBLEM_NAVIER_STOKES_EQUATION_TYPE

PROBLEM_SUBTYPE
PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE

END_PROBLEM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Creates a region to which geometric entities like ...!!!!!!!!!!!!!!!
!!!!!!!!! Cooridnate System, Mesh are attached !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_RegionType!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_REGION

REGION_ID
FLUID

REGION_LABEL_SET   ! THis is the group name that appreas in the .exelem ad .exnode file
Region

END_REGION



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Creates a Material Field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_FieldType!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_MATERIAL_FIELD

MATERIAL_FIELD_ID
FLUID

MATERIAL_FIELD_LABEL_SET  ! THis is the field name that appreas in the .exelem ad .exnode file
Material


MATERIAL_FIELD_PARAMETERS !! viscosity and density in this case
1,1                       !! user is supposed to know the order.


END_MATERIAL_FIELD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Creates a Geometric Field !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_FieldType!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_GEOMETRIC_FIELD

GEOMETRIC_FIELD_ID
FLUID

GEOMETRIC_FIELD_LABEL_SET  ! THis is the field name that appreas in the .exelem ad .exnode file
Geometry


END_GEOMETRIC_FIELD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines an incremental loop!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_ControlLoopType!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_CONTROL_LOOP

CONTROL_LOOP_ID
FLUID

CONTROL_LOOP_TYPE
PROBLEM_CONTROL_TIME_LOOP_TYPE

CONTROL_LOOP_TIME_INCREMENTS           !! initial_time,final_time,increment
0,0.1,0.001

! CONTROL_LOOP_LOAD_INCREMENTS         !!  THis option is used if the control loop type is
! 50                                   !!  PROBLEM_CONTROL_LOAD_LOOP_TYPE


END_CONTROL_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a solver !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_SolverType !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_SOLVER_SETTINGS

SOLVER_ID
FLUID

LINEAR_SOLVER_TYPE                     !! type of linear solver , direct or iterative
SOLVER_LINEAR_DIRECT_SOLVE_TYPE

SOLVER_LIBRARY_TYPE
SOLVER_MUMPS_LIBRARY

!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FEATURES/KEYWORDS SHOULD BE INCLUDED FOR THE DIRECT SOLVERS !!!!!!!!!!!!!!!!!!

EQUATION_DIRECT_SOLVER_TYPE
SOLVER_DIRECT_LU

!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FEATURES/KEYWORDS SHOULD BE INCLUDED FOR THE ITERATIVE SOLVERS !!!!!!!!!!!!!!!!!!

EQUATION_ITERATIVE_SOLVER_TYPE
SOLVER_ITERATIVE_GMRES

LINEAR_ITERATIVE_RELATIVE_TOLERANCE    !! User can can use both relative and absolute tolerance or ignore one of them
1E-6

LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE    !! User can can use both relative and absolute tolerance or ignore one of them
1E-12

LINEAR_ITERATIVE_MAXIMUM_ITERATIVE     !! Maximum iterations for the linear iterative solvers
100000

LINEAR_ITERATIVE_PRECONDITIONER_TYPE   !! THis parameter can be ignored by the user. The default type is no preconditioner
SOLVER_ITERATIVE_JACOBI_PRECONDITIONER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NONLINEAR_SOLVER                       !! for a quasi static problem this parameter is STATIC_NONLINEAR
Dynamic_nonlinear                      !! "_nonlinear" refers to the fact that the equilibrium is achieved for each load/time..
                                       !! .. increment using newton raphson iterations

NEWTON_SOLVER_MAXIMUM_ITERATIONS       !! Maximum number of Newton Raphson iterations
1000

NEWTON_SOLVER_ABSOLUTE_TOLERANCE       !! User can can use both relative and absolute tolerance or ignore one of them
1e-6

NEWTON_SOLVER_RELATIVE_TOLERANCE       !! User can can use both relative and absolute tolerance or ignore one of them
1e-6

JACOBIAN_TYPE                          !! default one is the one caluclated analytically using double derivation of the energy
SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED !! .. function

SOLVER_OUTPUT_TYPE_SET                 !! display solver output on screen m informations such as residual , total number of
SOLVER_PROGRESS_OUTPUT                 !! solver iterations

END_SOLVER_SETTINGS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a coordinate system !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_CoordinateSystemType!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_COORDINATE_SYSTEM

COORDINATE_SYSTEM_ID
FLUID

COORDINATE_SYSTEM_TYPE                 !! type cooridnate system
COORDINATE_RECTANGULAR_CARTESIAN_TYPE

COORDINATE_SYSTEM_DIMENSION            !! dimensions of the coordinate system
COORDINATE_SYSTEM_2D

END_COORDINATE_SYSTEM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a Dependent Field !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_FieldType!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_DEPENDENT_FIELD

DEPENDENT_FIELD_ID
FLUID


NUM_OF_COMPONENTS                      !! No. of state variables Vx, Vy and Pressure
3

INITIAL_VALUE_OF_VECTOR_STATE_VARIABLE !! Initial values for guess velocity for the first increment of Newton Raphson solver
0.0,0.0                                !! "undeformed" for an derformed eometry

INITIAL_VALUE_OF_SCALAR_STATE_VARIABLE !! Initial value of pressure state variable
0.0

DEPENDENT_FIELD_LABEL_SET              ! THis is the field name that appreas in the .exelem ad .exnode file
DependentField

END_DEPENDENT_FIELD




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a generated mesh !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_GeneratedMeshType!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_GENERATED_MESH

GENERATED_MESH_ID
FLUID


GENERATED_MESH_TYPE
GENERATED_REGULAR_MESH_TYPE

GENERATED_MESH_GEOMETRIC_EXTENTS    !! length of sides of rectangles
10,	2

GENERATED_MESH_ORIGIN_SET           !! location of left lower corner of the rectangle
0,0

GENERATED_MESH_NUMBER_OF_ELEMENTS   !!  no. of elements along each dimension
10,  10, 0                          !!  In a 2D case make sure to specify 0 elements in the third direction


END_GENERATED_MESH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a  mesh !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_MeshType !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_MESH

MESH_ID
FLUID

END_MESH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a velocity basis !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_BASIS !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_BASIS

BASIS_ID
FLUID

BASIS_INTERPOLATION_TYPE
QUADRATIC_LAGRANGE_INTERPOLATION

NumberOfGaussXi                                 !! no. of gauss point along each direction
3,3

END_BASIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a pressure basis !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_BASIS !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_BASIS


BASIS_ID
FLUID


NumberOfGaussXi
2  ,   2


BASIS_INTERPOLATION_TYPE
LINEAR_LAGRANGE_INTERPOLATION


END_BASIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines boundary condition!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_BoundaryCOnditionType!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_BOUNDARY_CONDITIONS

BOUNDARY_CONDITIONS_ID
FLUID

DIRICHELET  !!!! [ TYPE , LOCATION , COMPONENTS , VALUE ]

SURFACE,	MESH_REGULAR_LEFT_SURFACE,	100,		QUADRATIC_1  , BOUNDARY_CONDITION_FIXED
SURFACE,	MESH_REGULAR_LEFT_SURFACE,	010,		QUADRATIC_2  , BOUNDARY_CONDITION_FIXED
SURFACE,	MESH_REGULAR_FRONT_SURFACE,	110,		0            , BOUNDARY_CONDITION_FIXED
SURFACE,	MESH_REGULAR_BACK_SURFACE,	110,		0            , BOUNDARY_CONDITION_FIXED
SURFACE,	MESH_REGULAR_RIGHT_SURFACE,	001,		0            , BOUNDARY_CONDITION_FIXED



END_BOUNDARY_CONDITIONS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! FOLLOWING BLOCK DEFINES EXPORT PREFIX !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
START_OUTPUT

OUTPUT_ID
FLUID

NODE_EXPORT                     ! Name of .exnode file
Node_Data_fluid_flow

ELEMENT_EXPORT                  ! Name of exelem file
Element_Data_fluid_flow


END_OUTPUT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines boundary condition!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_DecompositionType!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
START_DECOMPOSITION

DECOMPOSITION_ID
FLUID

CALCULATE_ELEMENT_FACES
TRUE

END_DECOMPOSITION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! FIELDS BLOCK initilize object of a data type CMFE_FIELDS !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
START_FIELDS

FIELDS_ID
FLUID

END_FIELDS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINING Functions that are prescribed in BC block !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_FUNCTION

FUNCTION_ID
QUADRATIC_1

FUNCTION_CONSTANTS                               !!! A,B,C,E,F,G,H,I,J in equation ....
0,-4,0,0,0,0,0,8,0,0
                                                 !!! f(x,y,z) = A*x^2 + B*y^2 + C*z^2 + D*x*y + E*x*z + F*y*z ...
                                                 !!! ... + G*x + H*y + I*z + J
                                                 !!! The given constants correpond to a parabolic profile of ...
                                                 !! .. shape V = 4*V_max * (y/h)(1-(y/h)) , where V_max = 4 units and h is 2 units

END_FUNCTION

START_FUNCTION

FUNCTION_ID
QUADRATIC_2

FUNCTION_CONSTANTS                               !!! A,B,C,E,F,G,H,I,J in equation ....
0,0,0,0,0,0,0,0,0,0
                                                 !!! f(x,y,z) = A*x^2 + B*y^2 + C*z^2 + D*x*y + E*x*z + F*y*z ...
                                                 !!! ... + G*x + H*y + I*z + J
                                                 !!! The given constants correpond to zero BC


END_FUNCTION


STOP_PARSING

