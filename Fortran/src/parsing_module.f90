
!!!!!!!!!!!!!!!!!_______                     HEADER COMMENTS.             ___________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! THE MODULE ENCAPSULATED THE FOLLOWING INFORMATION: ..... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! _________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE parsing



  USE OpenCMISS
  USE OpenCMISS_Iron

  integer                                      :: i,FileStat,n,Stat,POS1,POS2,dummy                                             !! These variables are used as indices in the parsing algorithm...

  character(len=100)                           :: Rdline                                                                        !! This variable stores the lines READ  Statement.

  character(len=100),allocatable               :: EquationsSet_arg1(:,:),EquationsSet_arg2(:,:),EquationsSet_arg3(:,:), &       !! These data structures contain the input arguments given by...
                                                    & EquationsSet_arg4(:,:), EquationsSet_parsing(:)                           !! the user related to the EquationsSet.



  character(len=100),allocatable               :: Problem_arg1(:,:),Problem_arg2(:,:),Problem_arg3(:,:), &                      !! These data structures contain the input arguments given by...
                                                    & Problem_arg4(:,:), Problem_parsing(:)                                     !!  the user related to Problem.



  character(len=100),allocatable               :: MaterialField_arg1(:,:),MaterialField_arg2(:,:),MaterialField_arg3(:,:), &    !! These data structures contain the input arguments given by...
                                                    & MaterialField_parsing(:)                                                  !!  the user related to MaterialField.

  character(len=100),allocatable               :: Mesh_parsing(:),Mesh_arg1(:,:),Mesh_arg2(:,:),Mesh_arg3(:,:),Mesh_arg4(:,:), & !! These data structures contain the input arguments given by...
                                                    & Mesh_arg5(:,:)                                                            !! related to the mesh.

  character(len=100),allocatable               :: BC_parsing(:),BC_arg1(:,:),BC_arg2(:,:),BC_arg3(:,:),BC_arg4(:,:), &          !! related to the BoundaryCOndition.
                                                    & BC_arg5(:,:)


  character(len=100),allocatable               :: Solvers_arguments(:,:), Solvers_parsing(:)                                    !! related to the Solver.

  character(len=100),allocatable               :: Basis_arguments(:,:), Basis_parsing(:)                                        !! related to the BoundaryCOndition.


  character(len=100),allocatable               :: Region_arguments(:,:), Region_parsing(:)                                      !! related to the Region.

  character(len=100),allocatable               :: CoordinateSystem_arguments(:,:), CoordinateSystem_parsing(:)                  !! related to the CoordinateSystem.

  character(len=100),allocatable               :: ControlLoop_arguments(:,:), ControlLoop_parsing(:)                            !! related to the ControlLOOP.

  character(len=100),allocatable               :: DependentField_arguments(:,:), DependentField_parsing(:)                      !! related to the DepENDent FIeld.

  character(len=100),allocatable               :: FiberField_arg1(:,:), FiberField_arg2(:,:),FiberField_parsing(:)              !! related to the FiberFIeld.

  character(len=100),allocatable               :: PressureBasis_parsing(:), PressureBasis_arguments(:,:)                        !! related to the PressureBasis.

  character(len=100),allocatable               :: Output_parsing(:), Output_arguments(:,:)                                      !! related to the OUtput.

  character(len=100),allocatable               :: Field_parsing(:), Field_arg1(:,:),Field_arg2(:,:), Field_arg3(:,:)            !! related to the Field.

  character(len=100),allocatable               :: DecompositionParsing(:), DecompositionArguments(:,:)                        !! related to the Decompostion.
  contains

  !!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "FIELD  (GRAVITATIONAL FIELD etc. ) BLOCK" IN INPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!   NOTE THAT  THE REMAINING SUBROUTINES ALSO HAVE THE SAME STRUCTURE, THEREFORE .... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!    FOR EASE ONLY THE FOLLOWING SUBROUTINE IS COMMENTED. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   Field_parsing_subroutine(String1,Field_arg1,Field_arg2,Field_arg3,Rdline)


    integer                                            :: i,FileStat,n,Stat,k,POS1,POS2                              !! These variables are used as indices in the parsing algorithm.
    character(len=100) , dimension(:,:),allocatable    :: Field_arg1,Field_arg2,Field_arg3                           !! Arguments where the field related parameters are stored.
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1                                                    !! THis is the data structure that contains the keywords,
                                                                                                                     !! .... the algorithm looks for.



    IF (trim(Rdline)=="START_FIELD") then                                                                            !! Start of the Field block in the input file.

       Readline:DO

          READ(12,'(A)',IOSTAT=FileStat), Rdline
          Rdline =  strip_space(Rdline)                                                                              !! remove leading and triling spaces from the string.
          CALL remove_character(Rdline,"!")                                                                          !! remove comments off a string.

          IF (trim(Rdline)=="END_FIELD") EXIT                                                                        !! Terminate as soon as the cursor hits this keyword.

          IF (Rdline.EQ.String1(1)) then                                                                             !! Storing input argument 1 of Field block.
             READ(12,'(A)',IOSTAT=FileStat), Rdline
             CALL skip_blank_lines(Rdline)                                                                           !! This  ignores all the blank spaces.
             CALL remove_character(Rdline,"!")
             Rdline =  strip_space(Rdline)
             Field_arg1(:,1) = Rdline
          END IF


          IF (Rdline.EQ.String1(2)) then                                                                            !! Storing input argument 2 of Field block.
             READ(12,'(A)',IOSTAT=FileStat), Rdline
             CALL skip_blank_lines(Rdline)
             CALL remove_character(Rdline,"!")
             Rdline =  strip_space(Rdline)
             Field_arg2(:,1) = Rdline
          END IF

          IF (Rdline.EQ.String1(3)) then                                                                            !! Storing input argument 3 of Field block.
             READ(12,'(A)',IOSTAT=FileStat), Rdline
             CALL skip_blank_lines(Rdline)
             CALL remove_character(Rdline,"!")
             Rdline =  strip_space(Rdline)
             CALL split_inline_argument(Rdline,Field_arg3)                                                          !!Argument 3 is supposed to be a vector,therefore this subroutine is used to to store right

          END IF                                                                                                    !! values in the data structure  data.

      END DO Readline


    END IF                                                                                                          !! IF (trim(Rdline)=="START_FIELD")

  END subroutine Field_parsing_subroutine

  !!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "BOUDARY CONDITIONS"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   BC_parsing_subroutine(String1,boundary_conditions_arg1,boundary_conditions_arg2,boundary_conditions_arg3, &
                 & boundary_conditions_arg4, Rdline,num_of_dirichelet)


    integer                                            :: i,FileStat,n,Stat,k,POS1,POS2,num_of_dirichelet
    character(len=100) , dimension(:,:),allocatable    :: boundary_conditions_arg1,boundary_conditions_arg2, &
                                                            & boundary_conditions_arg4,boundary_conditions_arg3
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1



    IF (trim(Rdline)=="START_BOUNDARY_CONDITIONS") then
      num_of_dirichelet= 0
      num_of_traction_neumann=0
      num_of_pressure_neumann=0

      Readline: DO

        READ(12,'(A)'), Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")
        IF (trim(Rdline)=="END_BOUNDARY_CONDITIONS") EXIT

          IF (Rdline.EQ.String1(1)) then
             READ(12,'(A)',IOSTAT=FileStat), Rdline
             CALL skip_blank_lines(Rdline)
             CALL remove_character(Rdline,"!")
             Rdline =  strip_space(Rdline)
             boundary_conditions_arg1(:,1) = Rdline
          END IF

          IF (Rdline.EQ.String1(2)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            CALL skip_blank_lines(Rdline)
            CALL remove_character(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg2,num_of_dirichelet)
          ENDIF

          IF (Rdline.EQ.String1(3)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            CALL skip_blank_lines(Rdline)
            CALL remove_character(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg3,num_of_traction_neumann)
          END IF

          IF (Rdline.EQ.String1(4)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            CALL skip_blank_lines(Rdline)
            CALL remove_character(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg4,num_of_pressure_neumann)
          END IF


      END DO Readline

    END IF   !! IF (trim(Rdline)=="END_BOUNDARY_CONDITIONS")

  END subroutine BC_parsing_subroutine

  !!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "MATERIAL BLOCK"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   MaterialField_parsing_subroutine(String1,MaterialField_arg1,MaterialField_arg2,Rdline)


    integer                                            :: i,FileStat,n,Stat,k,POS1,POS2
    character(len=100) , dimension(:,:),allocatable    :: MaterialField_arg1,MaterialField_arg2
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1



    IF (trim(Rdline)=="START_MATERIAL_FIELD") then

      Readline:DO

        READ(12,'(A)',IOSTAT=FileStat), Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")

        IF (trim(Rdline)=="END_MATERIAL_FIELD") EXIT

        IF (Rdline.EQ.String1(1)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          MaterialField_arg1(:,1) = Rdline
        END IF


        IF (Rdline.EQ.String1(2)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          CALL split_inline_argument(Rdline,MaterialField_arg2)
        END IF

      END DO Readline

    END IF

  END subroutine MaterialField_parsing_subroutine



  !!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "FiberFIeld BLOCK"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   FiberField_parsing_subroutine(String1,FiberField_arg1,FiberField_arg2, Rdline)


    integer                                            :: i,FileStat,n,Stat,k,POS1,POS2
    character(len=100) , dimension(:,:),allocatable    :: FiberField_arg1,FiberField_arg2
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1



    IF (trim(Rdline)=="START_FIBER_FIELD") then


      Readline:DO

        READ(12,'(A)',IOSTAT=FileStat), Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")

        IF (trim(Rdline)=="END_FIBER_FIELD") EXIT

        IF (Rdline.EQ.String1(1)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          FiberField_arg1(:,1) = Rdline
        END IF

        IF (Rdline.EQ.String1(2)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          CALL split_inline_argument(Rdline,FiberField_arg2)
        END IF


      END DO Readline

    END IF

  END subroutine FiberField_parsing_subroutine

  !!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "Equation Set BLOCK"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   Equations_Set_parsing(String1,EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3, EquationsSet_arg4, Rdline)


    integer                                            :: i,FileStat,n,Stat,k
    character(len=100) , dimension(:,:),allocatable    :: EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3, &
                                                          & EquationsSet_arg4
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1



    IF (trim(Rdline)=="START_EQUATIONS_SET") then


      Readline:DO


        READ(12,'(A)',IOSTAT=FileStat), Rdline
        IF (trim(Rdline)=="END_EQUATIONS_SET") EXIT
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")

        IF (Rdline.EQ.String1(1)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          EquationsSet_arg1(:,1) = Rdline
        END IF

        IF (Rdline.EQ.String1(2)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          EquationsSet_arg2(:,1) = Rdline
        ENDIF

        IF (Rdline.EQ.String1(3)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          EquationsSet_arg3(:,1) = Rdline
        ENDIF

        IF (Rdline.EQ.String1(4)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          EquationsSet_arg4(:,1) = Rdline
        ENDIF


      END DO Readline

    END IF

  END subroutine Equations_Set_parsing

  !!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "PROBLEM BLOCK"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   Problem_parsing_subroutine(String1,Problem_arg1,Problem_arg2,Problem_arg3,Problem_arg4, Rdline)


    integer                                            :: i,FileStat,n,Stat,k
    character(len=100) , dimension(:,:),allocatable    :: Problem_arg1,Problem_arg2,Problem_arg3,Problem_arg4
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1




    IF (trim(Rdline)=="START_PROBLEM") then


      Readline:DO

        READ(12,'(A)',IOSTAT=FileStat), Rdline

        IF (trim(Rdline)=="END_PROBLEM") EXIT
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")

        IF (Rdline.EQ.String1(1)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          Rdline =  strip_space(Rdline)
          Problem_arg1(:,1) = Rdline
        END IF

        IF (Rdline.EQ.String1(2)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          Rdline =  strip_space(Rdline)
          Problem_arg2(:,1) = Rdline
        ENDIF

        IF (Rdline.EQ.String1(3)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          Rdline =  strip_space(Rdline)
          Problem_arg3(:,1) = Rdline
        ENDIF

        IF (Rdline.EQ.String1(4)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
                   Rdline =  strip_space(Rdline)
          Problem_arg4(:,1) = Rdline
        ENDIF

      END DO Readline

    END IF

  END subroutine Problem_parsing_subroutine

  !!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "MESH BLOCK"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine   Mesh_parsing_subroutine(String1,Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4,Rdline)

    integer                                            :: i,FileStat,n,Stat,k,POS1,POS2
    character(len=100) , dimension(:,:),allocatable    :: Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4 ,Mesh_arg5
    character(len=100)                                 :: Rdline,g
    character(len=100),   dimension(:)                 :: String1



    IF (trim(Rdline)=="START_MESH") then


      Readline:DO

        READ (12,'(A)',IOSTAT=FileStat), Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")


        IF (trim(Rdline)=="END_MESH") EXIT

        IF (Rdline.EQ.String1(1)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          Mesh_arg1(:,1) = Rdline
        END IF

        IF (Rdline.EQ.String1(2)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          Mesh_arg2(:,1) = Rdline
        END IF

        IF (Rdline.EQ.String1(3)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          CALL remove_character(Rdline,"!")
          Rdline =  strip_space(Rdline)
          CALL split_inline_argument(Rdline,Mesh_arg3)
        END IF


        IF (Rdline.EQ.String1(4)) then
          READ(12,'(A)',IOSTAT=FileStat), Rdline
          CALL skip_blank_lines(Rdline)
          Rdline =  strip_space(Rdline)
          CALL remove_character(Rdline,"!")
          k = 0
          POS1 = 1
          DO
            POS2 = INDEX(Rdline(POS1:), ",")
            IF (POS2 == 0) THEN
              k = k + 1
              Mesh_arg4(k,1) = Rdline(POS1:)
              EXIT
            END IF
            k = k + 1
            Mesh_arg4(k,1) = Rdline(POS1:POS1+POS2-2)
            POS1 = POS2+POS1
          END DO
        END IF

      END DO Readline

    END IF

  END subroutine Mesh_parsing_subroutine




  !!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE VERY BLOCK WITH SCALAR INPUTS     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE  parsing_subroutine(string,arguments,Rdline,num_of_arguments,string_tobe_parsed,ID)


    character(len=100),dimension(:,:)          :: arguments
    character(len=100)                         :: Rdline
    character(len=100),dimension(:)            :: string
    character(len=*)                           :: string_tobe_parsed
    integer                                    :: m,num_of_arguments



    IF (trim(Rdline)=="START_"// trim(string_tobe_parsed)) then
      ID = ID+1
      Readline:DO

        READ(12,'(A)',IOSTAT=FileStat), Rdline
        Rdline =  strip_space(Rdline)
        CALL remove_character(Rdline,"!")
        ! terminate upon encoutering END_start_string

        IF (trim(Rdline)=="END_" // string_tobe_parsed ) THEN
          EXIT
        END IF

        IF ( (trim(Rdline).NE."") .AND.  (trim(Rdline(:1)).NE."!")) then
          DO m = 1,num_of_arguments
            IF ((trim(Rdline)).EQ.string(m))  then
              READ(12,'(A)',IOSTAT=FileStat), arguments(m,ID)
              EXIT
            END IF
          END DO
        END IF

      END DO Readline
    END IF


  END SUBROUTINE  parsing_subroutine




  !!!!!!!!!!!!!!!!!!!!! THE FOLLOWING SUBROUTINE COUNTS THE NUMBER OF TIMES A KEYWORD EXIST IN THE INPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE  searching(Rdline,search_string,counter)

    !!!! initializing values


    character(len=100)                         :: Rdline
    character(len=*)                           :: search_string
    integer                                    :: counter

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
    INTEGER :: iachar_character

    ! -- Initialise output string
    output_string = ' '

    ! -- Initialise output string "useful" length counter
    n = 0

    ! -- Loop over string elements
    DO i = 1, LEN( input_string )

      ! -- Convert the current character to its position
      ! -- in the ASCII collating sequence
      iachar_character = IACHAR( input_string( i:i ) )

      ! -- IF the character is NOT a space ' ' or a tab '->|'
      ! -- copy it to the output string.
      IF ( iachar_character /= IACHAR_SPACE .AND. &
           iachar_character /= IACHAR_TAB) THEN
        n = n + 1
        output_string( n:n ) = input_string( i:i )
      END IF

    END DO


  END FUNCTION strip_space



  !!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION... !!!!!!!!!!!!!!!!!!!!!!
  !!!!!         REMOVE SPECIFIED CHARACTER FROM A STRING!!!!!!!!!!!!!!!!



  subroutine remove_character(string,character)
    implicit none
    integer :: ind
    character(len=100)  :: string
    character           :: character
    ind = SCAN(string, character)
    IF (ind .NE.0) then
      string = string(:ind-1)
    END IF
  END subroutine remove_character




  !!!!!CONVERT STRING TO REAL NUMBER !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!



  REAL(CMISSRP) function str2real(str_2)
    implicit none
    ! Arguments
    character(len=*)      :: str_2
    READ(str_2,*,IOSTAT=Stat)  str2real
  END function str2real


  !!!!!CONVERT STRING TO INTEGER NUMBER !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!

  integer function str2int(str)
    implicit none
    ! Arguments
    character(len=*)      :: str
    READ(str,*,IOSTAT=Stat)  str2int
  END function str2int

  !!!!!Skip blank lines in the input file !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!

  subroutine skip_blank_lines(Rdline)


    integer :: FileStat
    character(len=100)              :: Rdline

    DO
      IF (Rdline .NE."") EXIT
      READ(12,*) Rdline
    END DO

  END subroutine skip_blank_lines


  !!!!!!!!!!!!!!!!!!!!!!!!!! The following subroutine splits a comma separated string !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine split_inline_argument(Rdline,array1,i)


    integer :: k,POS2,POS1
    integer,optional::i
    character(len=100)  :: Rdline
    character(len=100) , dimension(:,:),allocatable :: array1
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
        CALL remove_character(Rdline,"!")
        IF(present(i)) i=i+1
   END DO

  END subroutine split_inline_argument


END MODULE parsing
