
!!!!!!!!!!!!!!!!!_______                     HEADER COMMENTS.             ___________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! THE MODULE ENCAPSULATED THE FOLLOWING INFORMATION: ..... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! _________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE parsing



  USE OpenCMISS
  USE OpenCMISS_Iron

  INTEGER                                      :: i,FileStat,n,Stat,POS1,POS2,dummy                                             !! These variables are used as indices in the parsing algorithm...

  CHARACTER(LEN=100)                           :: Rdline                                                                        !! This variable stores the lines READ  Statement.

  CHARACTER(LEN=100),ALLOCATABLE               :: EquationSetKeyWords(:),ProblemKeyWords(:),MaterialFieldKeyWords(:)            !! Data Structures holding keywords that the parsing algorithm...
  CHARACTER(LEN=100),ALLOCATABLE               :: FibreFieldKeyWords(:),SolverKeywords(:),CoordinateSystemKeywords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: MeshKeyWords(:),DecompositionKeyWords(:),DependentFieldKeyWords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: BasisKeyWords(:),OutputKeyWords(:),BoundaryConditionKeywords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: RegionKeywords(:),ControlLoopKeyWords(:),FieldKeyWords(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES TO STORE USER INPUTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TYPE ParsinDataStructures

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: EquationSetId(:,:),EquationSetClass(:,:),EquationSetType(:,:), &
                                                    & EquationSetSubType(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: ProblemId(:,:),ProblemClass(:,:),ProblemType(:,:),ProblemSubType(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: MaterialFieldId(:,:),MaterialFieldParameters(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: FibreFieldId(:,:),FibreFieldParameters(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: SolverId(:,:),SolverType(:,:),Preconditioner(:,:), &
                                                    & EquationSolverMaximumIteration(:,:), &
                                                      & EquationSolverTolerance(:,:),NewtonMaxIterations(:,:), &
                                                        & NewtonTolerance(:,:), NewtonJacobianType(:,:), &
                                                          & NonlinearSolver(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: CoordinateSystemId(:,:),CoordinateSystemType(:,:),CoordinateSystemDimension(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: MeshID(:,:),MeshTopology(:,:),MeshGeometricParameters(:,:), &
                                                    & MeshNumberOfElements(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: NumberOfDomains(:,:),DecompositionFaceActive(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: DependentFieldId(:,:),DependentFieldStateVariable(:,:), &
                                                    & DependentFieldNumberOfComponents(:,:), &
                                                      & DependentFieldInitialValueOfStateVector(:,:), &
                                                        & DependentFieldInitialValueOfStateScalar(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: BasisId(:,:),BasisType(:,:),BasisComponents(:,:), &
                                                    & BasisNumberOfGaussPoints(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: ControlLoopId(:,:),ControlLoopType(:,:),ControlLoopTimeIncrement(:,:), &
                                                    & ControlLoopLoadIncrement(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: MatrixType(:,:),OutputType(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: RegionId(:,:)

    TYPE(CHARACTER(LEN=100)), ALLOCATABLE      :: SourceFieldId(:,:) , SourceFieldType(:,:), &
                                                  SourceFieldComponents(:,:)



  END TYPE ParsinDataStructures

    TYPE(ParsinDataStructures) 	               :: DependentField,Problem,MaterialField,Mesh
    TYPE(ParsinDataStructures) 	               :: BoundaryCondition,CoordinateSystem,Region
    TYPE(ParsinDataStructures) 	               :: Solver,Basis,ControlLoop,FibreField,SourceField
    TYPE(ParsinDataStructures) 	               :: Output,Field,Decomposition,EquationsSet
    CHARACTER(LEN=100),ALLOCATABLE             :: BoundaryConditionId(:,:),BoundaryDirichelet(:,:), &
                                                    & BoundaryTractionNeumann(:,:), BoundaryPressureNeumann(:,:)
  CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!                              EXPERIMENT BLOCK                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  SUBROUTINE   GenericParsing(KeywordArray, BlockToBeParsed,NumOfArguments,Id,Rdline, Argument1,Argument2,Argument3, &
                                & Argument4,Argument5,Argument6,Argument7,Argument8,Argument9)


    INTEGER                                                    :: counter,FileStat,n,Stat,k,POS1,POS2,NumOfArguments,Id       !! These variables are used as indices in the parsing algorithm.
    CHARACTER(LEN=*)    , OPTIONAL   , DIMENSION(:)            :: Argument1,Argument2,Argument3,Argument4,Argument5
    CHARACTER(LEN=*)    , OPTIONAL   , DIMENSION(:)            :: Argument6,Argument7,Argument8,Argument9                     !! Arguments where the field related parameters are stored.
    CHARACTER(LEN=*)                                           :: Rdline
    CHARACTER(LEN=*)    , DIMENSION(:)                         :: KeywordArray                                                !! THis is the data structure that contains the keywords,
                                                                                                                              !! .... the algorithm looks for within a block.
    CHARACTER(LEN=*)                                           :: BlockToBeParsed                                             !! Block to be parsed in the input file.
    CHARACTER(LEN=100)  , ALLOCATABLE                          :: StoreArguments(:,:)
    LOGICAL                                                    :: Flag

    ALLOCATE(StoreArguments(9,9))

    IF (TRIM(Rdline)== "START_" // TRIM(BlockToBeParsed)) THEN        !! Start of the Field block in the input file.

       Id = Id + 1
       Readline:DO

          counter=0
          Flag = .TRUE.
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          Rdline =  strip_space(Rdline)
          CALL remove_CHARACTER(Rdline,"!")
          IF (TRIM(Rdline)=="END_"//TRIM(BlockToBeParsed)) THEN
            EXIT
          END IF


          DO WHILE (Flag .EQV. .TRUE. .AND. counter .LE. NumOfArguments)
            counter = counter + 1

            IF (Rdline.EQ.KeywordArray(counter)) THEN
              Flag = .FALSE.
              READ(12,'(A)',IOSTAT=FileStat), Rdline
              CALL skip_blank_lines(Rdline)
              CALL remove_CHARACTER(Rdline,"!")
              Rdline =  strip_space(Rdline)

              CALL split_inline_argument_new(Rdline,StoreArguments(counter,:))

            END IF
          END DO

       !! values in the data structure  data.
      END DO Readline

       Argument1 = StoreArguments(1,:)


       IF  (NumOfArguments .GT. 1 ) THEN
         Argument2 = StoreArguments(2,:)
       END IF

       IF  (NumOfArguments .GT. 2 ) THEN
         Argument3 = StoreArguments(3,:)
       END IF

       IF  (NumOfArguments .GT. 3 ) THEN
         Argument4 = StoreArguments(4,:)

       END IF

       IF  (NumOfArguments .GT. 4 ) THEN
         Argument5 = StoreArguments(5,:)
       END IF

       IF  (NumOfArguments .GT. 5 ) THEN
         Argument6 = StoreArguments(6,:)
       END IF

       IF  (NumOfArguments .GT. 6 ) THEN
         Argument7 = StoreArguments(7,:)
       END IF

       IF  (NumOfArguments .GT. 7 ) THEN
         Argument8 = StoreArguments(8,:)
       END IF

       IF  (NumOfArguments .GT. 8 ) THEN
         Argument9 = StoreArguments(9,:)
       END IF

    END IF


  END SUBROUTINE GenericParsing


  !!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "BOUDARY CONDITIONS"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE   BC_parsing_subroutine(String1,boundary_conditions_arg1,boundary_conditions_arg2,boundary_conditions_arg3, &
                 & boundary_conditions_arg4, Rdline,num_of_dirichelet)


    INTEGER                                            :: i,FileStat,n,Stat,k,POS1,POS2,num_of_dirichelet
    CHARACTER(LEN=100) , dimension(:,:),ALLOCATABLE    :: boundary_conditions_arg1,boundary_conditions_arg2, &
                                                            & boundary_conditions_arg4,boundary_conditions_arg3
    CHARACTER(LEN=100)                                 :: Rdline,g
    CHARACTER(LEN=100),   dimension(:)                 :: String1



    IF (trim(Rdline)=="START_BOUNDARY_CONDITIONS") then
      num_of_dirichelet= 0
      num_of_traction_neumann=0
      num_of_pressure_neumann=0

      Readline: DO

        READ(12,'(A)'), Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_CHARACTER(Rdline,"!")
        IF (trim(Rdline)=="END_BOUNDARY_CONDITIONS") EXIT

          IF (Rdline.EQ.String1(1)) then
             READ(12,'(A)',IOSTAT=FileStat), Rdline
             CALL skip_blank_lines(Rdline)
             CALL remove_CHARACTER(Rdline,"!")
             Rdline =  strip_space(Rdline)
             boundary_conditions_arg1(:,1) = Rdline
          END IF

          IF (Rdline.EQ.String1(2)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            CALL skip_blank_lines(Rdline)
            CALL remove_CHARACTER(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg2,num_of_dirichelet)
          ENDIF

          IF (Rdline.EQ.String1(3)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            CALL skip_blank_lines(Rdline)
            CALL remove_CHARACTER(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg3,num_of_traction_neumann)
          END IF

          IF (Rdline.EQ.String1(4)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            CALL skip_blank_lines(Rdline)
            CALL remove_CHARACTER(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg4,num_of_pressure_neumann)
          END IF


      END DO Readline

    END IF   !! IF (trim(Rdline)=="END_BOUNDARY_CONDITIONS")

  END SUBROUTINE BC_parsing_subroutine

  !!!!!!!!!!!!!!!!!!!!! THE FOLLOWING SUBROUTINE COUNTS THE NUMBER OF TIMES A KEYWORD EXIST IN THE INPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE  searching(Rdline,search_string,counter)

    !!!! initializing values


    CHARACTER(LEN=100)                         :: Rdline
    CHARACTER(LEN=*)                           :: search_string
    INTEGER                                    :: counter

    Rdline =  strip_space(Rdline)
    IF (trim(Rdline)==search_string) counter = counter+1

  END SUBROUTINE  searching




  !!!!!!!!!!!REMOVE PROCEEDING , RECEEDING SPACES AND SPACES INBETWEEN CHARACTERS INSIDE A STRING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION strip_space( input_string ) RESULT ( output_string )

    ! -- Arguments
    CHARACTER( * ), INTENT( IN )  :: input_string
    INTEGER                       :: n

    ! -- Function result
    CHARACTER( LEN( input_string ) ) :: output_string

    ! -- Local parameters
    INTEGER,        PARAMETER :: IACHAR_SPACE = 32, IACHAR_TAB   = 9


    ! -- Local variables
    INTEGER :: i
    INTEGER :: iachar_CHARACTER

    ! -- Initialise output string
    output_string = ' '

    ! -- Initialise output string "useful" LENgth counter
    n = 0

    ! -- Loop over string elements
    DO i = 1, LEN( input_string )

      ! -- Convert the current CHARACTER to its position
      ! -- in the ASCII collating sequence
      iachar_CHARACTER = IACHAR( input_string( i:i ) )

      ! -- IF the CHARACTER is NOT a space ' ' or a tab '->|'
      ! -- copy it to the output string.
      IF ( iachar_CHARACTER /= IACHAR_SPACE .AND. &
           iachar_CHARACTER /= IACHAR_TAB) THEN
        n = n + 1
        output_string( n:n ) = input_string( i:i )
      END IF

    END DO


  END FUNCTION strip_space



  !!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION... !!!!!!!!!!!!!!!!!!!!!!
  !!!!!         REMOVE SPECIFIED CHARACTER FROM A STRING!!!!!!!!!!!!!!!!



  SUBROUTINE remove_CHARACTER(string,CHARACTER)
    implicit none
    INTEGER :: ind
    CHARACTER(LEN=100)  :: string
    CHARACTER           :: CHARACTER
    ind = SCAN(string, CHARACTER)
    IF (ind .NE.0) then
      string = string(:ind-1)
    END IF
  END SUBROUTINE remove_CHARACTER




  !!!!!CONVERT STRING TO REAL NUMBER !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!



  REAL(CMISSRP) function str2real(str_2)
    implicit none
    ! Arguments
    CHARACTER(LEN=*)      :: str_2
    READ(str_2,*,IOSTAT=Stat)  str2real
  END function str2real


  !!!!!CONVERT STRING TO INTEGER NUMBER !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!

  INTEGER function str2int(str)
    implicit none
    ! Arguments
    CHARACTER(LEN=*)      :: str
    READ(str,*,IOSTAT=Stat)  str2int
  END function str2int

  !!!!!Skip blank lines in the input file !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!

  SUBROUTINE skip_blank_lines(Rdline)


    INTEGER :: FileStat
    CHARACTER(LEN=100)              :: Rdline

    DO
      IF (Rdline .NE."") EXIT
      READ(12,*) Rdline
    END DO

  END SUBROUTINE skip_blank_lines

    !!!!!!!!!!!!!!!!!!!!!!!!!! The following SUBROUTINE splits a comma separated string !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE split_inline_argument_new(Rdline,array1,i)


    INTEGER :: k,POS2,POS1
    INTEGER,optional::i
    CHARACTER(LEN=100)  :: Rdline
    CHARACTER(LEN=*) , dimension(:) :: array1
    k=0
    DO

        POS1 = 1
        DO
            POS2 = INDEX(Rdline(POS1:), ",")
            IF (POS2 == 0) THEN
                k = k + 1
                array1(k) = Rdline(POS1:)
                EXIT
            END IF
            k = k + 1
            array1(k) = Rdline(POS1:POS1+POS2-2)
            POS1 = POS2+POS1
        END DO

        IF (Rdline == "") EXIT
        READ(12,'(A)') Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_CHARACTER(Rdline,"!")
        IF(present(i)) i=i+1

   END DO

  END SUBROUTINE split_inline_argument_new


  !!!!!!!!!!!!!!!!!!!!!!!!!! The following SUBROUTINE splits a comma separated string !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE split_inline_argument(Rdline,array1,i)


    INTEGER :: k,POS2,POS1
    INTEGER,optional::i
    CHARACTER(LEN=100)  :: Rdline
    CHARACTER(LEN=100) , dimension(:,:),ALLOCATABLE :: array1
    k=0
    DO

        POS1 = 1
        DO
            POS2 = INDEX(Rdline(POS1:), ",")
            IF (POS2 == 0) THEN
                k = k + 1
                array1(k,1) = Rdline(POS1:)
                EXIT
            END IF
            k = k + 1
            array1(k,1) = Rdline(POS1:POS1+POS2-2)
            POS1 = POS2+POS1
        END DO

        IF (Rdline == "") EXIT
        READ(12,'(A)') Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_CHARACTER(Rdline,"!")
        IF(present(i)) i=i+1
   END DO

  END SUBROUTINE split_inline_argument


END MODULE parsing
