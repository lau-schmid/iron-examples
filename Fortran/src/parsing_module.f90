module parsing

  USE OpenCMISS
  USE OpenCMISS_Iron

  CHARACTER(*), PARAMETER                      :: fileplace = "/data/homes/mirzawd/software/opencmiss_examples/&
                                                              &large_uniaxial_extension/Fortran/src/"

  integer                                      :: i,filestat,n,stat,pos1,pos2
  character(len=100)                           :: rdline
  character(len=100),allocatable               :: EquationsSet_arg1(:,:),EquationsSet_arg2(:,:),EquationsSet_arg3(:,:), &
                                                &EquationsSet_arg4(:,:), EquationsSet_parsing(:)

  character(len=100),allocatable               :: Problem_arg1(:,:),Problem_arg2(:,:),Problem_arg3(:,:), &
                                                &Problem_arg4(:,:), Problem_parsing(:)


  character(len=100),allocatable               :: MaterialField_arg1(:,:),MaterialField_arg2(:,:),MaterialField_arg3(:,:), &
                                                &MaterialField_parsing(:)

  character(len=100),allocatable               :: Mesh_parsing(:),Mesh_arg1(:,:),Mesh_arg2(:,:),Mesh_arg3(:,:),Mesh_arg4(:,:),&
                                                 & Mesh_arg5(:,:)

  character(len=100),allocatable               :: BC_parsing(:),BC_arg1(:,:),BC_arg2(:,:),BC_arg3(:,:),BC_arg4(:,:),&
                                                &BC_arg5(:,:)


  character(len=100),allocatable               :: Solvers_arguments(:,:), Solvers_parsing(:)

  character(len=100),allocatable               :: Basis_arguments(:,:), Basis_parsing(:)


   character(len=100),allocatable              :: Region_arguments(:,:), Region_parsing(:)

   character(len=100),allocatable              :: CoordinateSystem_arguments(:,:), CoordinateSystem_parsing(:)

   character(len=100),allocatable              :: ControlLoop_arguments(:,:), ControlLoop_parsing(:)

   character(len=100),allocatable              :: DependentField_arguments(:,:), DependentField_parsing(:)

   character(len=100),allocatable              :: FiberField_arg1(:,:), FiberField_arg2(:,:),FiberField_parsing(:)

contains


 subroutine   BC_parsing_subroutine(string1,boundary_conditions_arg1,boundary_conditions_arg2,boundary_conditions_arg3,&
                       & boundary_conditions_arg4, rdline)

 integer                                            :: i,filestat,n,stat,k,pos1,pos2
 character(len=100) , dimension(:,:),allocatable    :: boundary_conditions_arg1,boundary_conditions_arg2,&
                                                      &boundary_conditions_arg4,boundary_conditions_arg3
 character(len=100)                                 :: rdline,g
 character(len=100),   dimension(:)                 :: string1

 k=0


 if (trim(rdline)=="START_BOUNDARY_CONDITIONS") then


 readline:do

    read(12,'(A)',iostat=filestat), rdline
   ! ! print *, 456
    !!!! TERMINATE UPON ENDIF

    if (trim(rdline)=="END_BOUNDARY_CONDITIONS") then
 !    ! print *, "i am here "
 !    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg3
!    ! print *, boundary_conditions_arg4
    exit
    end if


    rdline =  strip_space(rdline)
    call remove_character(rdline,"!")

    if ( (trim(rdline).NE."") .AND. (rdline(:1).NE."!")) then


    if (rdline.EQ.string1(1)) then

        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        boundary_conditions_arg1(:,1) = rdline

 !       ! print *, boundary_conditions_arg1(1,1)

    end if

    if (rdline.EQ.string1(2)) then
 do i = 1,4
        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")

!       ! print *, rdline
!        call split_string(rdline,boundary_conditions_arg2(:,1))
        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       boundary_conditions_arg2(k+4*(i-1),1) = rdline(pos1:)
    !   print  *,  boundary_conditions_arg2(k*i,1)
       EXIT
    END IF
    k = k + 1
    boundary_conditions_arg2(k+4*(i-1),1) = rdline(pos1:pos1+pos2-2)
  !  print  *,  boundary_conditions_arg2(k*i,1)
    pos1 = pos2+pos1

   END DO

end do

   endif
! ! print *, rdline



    if (rdline.EQ.string1(3)) then
 do i = 1,1
        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")

!       ! print *, rdline
!        call split_string(rdline,boundary_conditions_arg2(:,1))
        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       boundary_conditions_arg3(k+4*(i-1),1) = rdline(pos1:)
  !     print  *, boundary_conditions_arg3(k*i,1)
       EXIT
    END IF
    k = k + 1
    boundary_conditions_arg3(k+4*(i-1),1) = rdline(pos1:pos1+pos2-2)
    pos1 = pos2+pos1

   END DO

end do

end if

  if (rdline.EQ.string1(4)) then
 do i = 1,1
        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")

!       ! print *, rdline
!        call split_string(rdline,boundary_conditions_arg2(:,1))
        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       boundary_conditions_arg4(k+4*(i-1),1) = rdline(pos1:)
!       print  *, boundary_conditions_arg4(k,1)
       EXIT
    END IF
    k = k + 1
    boundary_conditions_arg4(k+4*(i-1),1) = rdline(pos1:pos1+pos2-2)
    pos1 = pos2+pos1

   END DO

end do


    end if



 end if


end do readline

end if
end subroutine BC_parsing_subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine   MaterialField_parsing_subroutine(string1,MaterialField_arg1,MaterialField_arg2,&
                                              & rdline)

 integer                                            :: i,filestat,n,stat,k,pos1,pos2
 character(len=100) , dimension(:,:),allocatable    :: MaterialField_arg1,MaterialField_arg2
 character(len=100)                                 :: rdline,g
 character(len=100),   dimension(:)                 :: string1
 k=0


 if (trim(rdline)=="START_MATERIAL_FIELD") then


 readline:do

    read(12,'(A)',iostat=filestat), rdline
    rdline =  strip_space(rdline)

    call remove_character(rdline,"!")


    if (trim(rdline)=="END_MATERIAL_FIELD") then
 !    ! print *, "i am here "
 !    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg3
!    ! print *, boundary_conditions_arg4
    exit
    end if


    rdline =  strip_space(rdline)
    call remove_character(rdline,"!")

    if ( (trim(rdline).NE."") .AND. (rdline(:1).NE."!")) then


    if (rdline.EQ.string1(1)) then

        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")
        MaterialField_arg1(:,1) = rdline
 !       ! print *, MaterialField_arg1(:,1)
 !       ! print *, boundary_conditions_arg1(1,1)

    end if


! ! print *, rdline



    if (rdline.EQ.string1(2)) then



        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")


        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       MaterialField_arg2(k,1) = rdline(pos1:)
! ! print *,  MaterialField_arg2(k,1)
       EXIT
    END IF
    k = k + 1
    MaterialField_arg2(k,1) = rdline(pos1:pos1+pos2-2)
    pos1 = pos2+pos1
! ! print *, MaterialField_arg2(k,1)
   END DO

end if

end if


end do readline

end if

end subroutine MaterialField_parsing_subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine   FieldField_parsing_subroutine(string1,FiberField_arg1,FiberField_arg2,&
                                              & rdline)

 integer                                            :: i,filestat,n,stat,k,pos1,pos2
 character(len=100) , dimension(:,:),allocatable    :: FiberField_arg1,FiberField_arg2
 character(len=100)                                 :: rdline,g
 character(len=100),   dimension(:)                 :: string1
 k=0


 if (trim(rdline)=="START_FIBER_FIELD") then


 readline:do

    read(12,'(A)',iostat=filestat), rdline
    rdline =  strip_space(rdline)

    call remove_character(rdline,"!")


    if (trim(rdline)=="END_FIBER_FIELD") then
 !    ! print *, "i am here "
 !    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg3
!    ! print *, boundary_conditions_arg4
    exit
    end if


    rdline =  strip_space(rdline)
    call remove_character(rdline,"!")

    if ( (trim(rdline).NE."") .AND. (rdline(:1).NE."!")) then


    if (rdline.EQ.string1(1)) then

        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")
        FiberField_arg1(:,1) = rdline
 !       ! print *, FiberField_arg1(:,1)
 !       ! print *, boundary_conditions_arg1(1,1)

    end if


! ! print *, rdline



    if (rdline.EQ.string1(2)) then



        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")


        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       FiberField_arg2(k,1) = rdline(pos1:)
! ! print *,  FiberField_arg2(k,1)
       EXIT
    END IF
    k = k + 1
    FiberField_arg2(k,1) = rdline(pos1:pos1+pos2-2)
    pos1 = pos2+pos1
! print *, FiberField_arg2(k,1)
   END DO

end if

end if


end do readline

end if

end subroutine FieldField_parsing_subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine   Equations_Set_parsing(string1,EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3, &
                                 &EquationsSet_arg4, rdline)

 integer                                            :: i,filestat,n,stat,k
 character(len=100) , dimension(:,:),allocatable    :: EquationsSet_arg1,EquationsSet_arg2,EquationsSet_arg3,&
                                                      &EquationsSet_arg4
 character(len=100)                                 :: rdline,g
 character(len=100),   dimension(:)                 :: string1
 k=0


 if (trim(rdline)=="START_EQUATIONS_SET") then


 readline:do


    read(12,'(A)',iostat=filestat), rdline

    !!!! TERMINATE UPON ENDIF
! ! print *, trim(rdline)=="START_EQUATIONS_SET"


    if (trim(rdline)=="END_EQUATIONS_SET") then
      ! print *, "i am here "
      ! print *, EquationsSet_arg1(:,1)
      ! print *, EquationsSet_arg2(:,1)
      ! print *, EquationsSet_arg3(:,1)
      ! print *, EquationsSet_arg4(:,1)
        exit

    end if


    rdline =  strip_space(rdline)
    call remove_character(rdline,"!")

    if ( (trim(rdline).NE."") .AND. (rdline(:1).NE."!")) then


    if (rdline.EQ.string1(1)) then


        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        EquationsSet_arg1(:,1) = rdline

 !       ! print *, boundary_conditions_arg1(1,1)

    end if

    if (rdline.EQ.string1(2)) then

        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        EquationsSet_arg2(:,1) = rdline


    endif
! ! print *, rdline

    if (rdline.EQ.string1(3)) then

        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        EquationsSet_arg3(:,1) = rdline


    endif



    if (rdline.EQ.string1(4)) then

        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        EquationsSet_arg4(:,1) = rdline


    endif

end if
end do readline
end if

end subroutine Equations_Set_parsing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine   Problem_parsing_subroutine(string1,Problem_arg1,Problem_arg2,Problem_arg3, &
                                       &Problem_arg4, rdline)

 integer                                            :: i,filestat,n,stat,k

 character(len=100) , dimension(:,:),allocatable    :: Problem_arg1,Problem_arg2,Problem_arg3,Problem_arg4

 character(len=100)                                 :: rdline,g

 character(len=100),   dimension(:)                 :: string1

 k=0


 if (trim(rdline)=="START_PROBLEM") then


 readline:do


    read(12,'(A)',iostat=filestat), rdline

    !!!! TERMINATE UPON ENDIF
! ! print *, trim(rdline)=="START_EQUATIONS_SET"


    if (trim(rdline)=="END_PROBLEM") then
      ! print *, Problem_arg1(:,1)
      ! print *, Problem_arg2(:,1)
      ! print *, Problem_arg3(:,1)
      ! print *, Problem_arg4(:,1)

        exit
    end if


    rdline =  strip_space(rdline)
    call remove_character(rdline,"!")

    if ( (trim(rdline).NE."") .AND. (rdline(:1).NE."!")) then


    if (rdline.EQ.string1(1)) then


        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        Problem_arg1(:,1) = rdline

 !       ! print *, boundary_conditions_arg1(1,1)

    end if

    if (rdline.EQ.string1(2)) then

        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        Problem_arg2(:,1) = rdline


    endif
! ! print *, rdline

    if (rdline.EQ.string1(3)) then

        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        Problem_arg3(:,1) = rdline


    endif



    if (rdline.EQ.string1(4)) then

        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        Problem_arg4(:,1) = rdline


    endif

end if
end do readline
end if

end subroutine Problem_parsing_subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine   Mesh_parsing_subroutine(string1,Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4,rdline)

 integer                                            :: i,filestat,n,stat,k,pos1,pos2
 character(len=100) , dimension(:,:),allocatable    :: Mesh_arg1,Mesh_arg2,Mesh_arg3,Mesh_arg4 ,Mesh_arg5
 character(len=100)                                 :: rdline,g
 character(len=100),   dimension(:)                 :: string1
 k=0


 if (trim(rdline)=="START_MESH") then


 readline:do

    read(12,'(A)',iostat=filestat), rdline
    rdline =  strip_space(rdline)

    call remove_character(rdline,"!")


    if (trim(rdline)=="END_MESH") then
 !    ! print *, "i am here "
 !    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg2
!    ! print *, boundary_conditions_arg3
!    ! print *, boundary_conditions_arg4
    exit
    end if


    rdline =  strip_space(rdline)
    call remove_character(rdline,"!")

    if ( (trim(rdline).NE."") .AND. (rdline(:1).NE."!")) then


    if (rdline.EQ.string1(1)) then

        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")
        Mesh_arg1(:,1) = rdline
        ! print *, Mesh_arg1(:,1)
 !       ! print *, boundary_conditions_arg1(1,1)

    end if

    if (rdline.EQ.string1(2)) then

        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline
        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")
        Mesh_arg2(:,1) = rdline
        ! print *, Mesh_arg2(:,1)
 !       ! print *, boundary_conditions_arg1(1,1)

    end if
! ! print *, rdline



    if (rdline.EQ.string1(3)) then



        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")


        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       Mesh_arg3(k,1) = rdline(pos1:)

       EXIT
    END IF
    k = k + 1
    Mesh_arg3(k,1) = rdline(pos1:pos1+pos2-2)
    pos1 = pos2+pos1

   END DO

end if


     if (rdline.EQ.string1(4)) then



        ! ! print *, rdline
        read(12,'(A)',iostat=filestat), rdline

        rdline =  strip_space(rdline)

        call remove_character(rdline,"!")


        k = 0
        pos1 = 1

    DO
    pos2 = INDEX(rdline(pos1:), ",")

    IF (pos2 == 0) THEN

       k = k + 1

       Mesh_arg4(k,1) = rdline(pos1:)

       EXIT
    END IF
    k = k + 1
    Mesh_arg4(k,1) = rdline(pos1:pos1+pos2-2)
    pos1 = pos2+pos1

   END DO

end if

end if


end do readline

end if
end subroutine Mesh_parsing_subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  Solver_parsing_subroutine(string,arguments,rdline)

!!!! initializing values
 Implicit None

 character(len=100)                         :: arguments(:,:)
 character(len=100)                         :: rdline
 character(len=100)                         :: string(:)
 integer                                    :: m

 if (trim(rdline)=="START_SOLVER_SETTINGS") then

 readline:do

        read(12,'(A)',iostat=filestat), rdline
        ! terminate upon encoutering END_start_string
        if (trim(rdline)=="END_SOLVER_SETTINGS" ) then

                exit

        end if

! to remove spaces
        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")
 !!!! TO IGNORE BLANK SPACES AND COMMENTS

        if ( (trim(rdline).NE."") .AND.  (trim(rdline(:1)).NE."!")) then
!!! search for key words and store values


                  do m = 1,8



                         if ((trim(rdline)).EQ.string(m))  then
                                read(12,'(A)',iostat=filestat), arguments(m,1)
                                ! print *, arguments(m,1)
                                exit
                        end if



                end do

         end if


 !!! MISCELLANEOUS STUFF

 end do readline
 end if


END SUBROUTINE  solver_parsing_subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!!!!!!!
SUBROUTINE  Basis_parsing_subroutine(string,arguments,rdline)

!!!! initializing values
 Implicit None

 character(len=100)                         :: arguments(:,:)
 character(len=100)                         :: rdline
 character(len=100)                         :: string(:)
 integer                                    :: m

 if (trim(rdline)=="START_BASIS") then

 readline:do

        read(12,'(A)',iostat=filestat), rdline
        ! terminate upon encoutering END_start_string
        if (trim(rdline)=="END_BASIS" ) then

                exit

        end if

! to remove spaces

        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")

 !!!! TO IGNORE BLANK SPACES AND COMMENTS

        if ( (trim(rdline).NE."") .AND.  (trim(rdline(:1)).NE."!")) then
!!! search for key words and store values


                  do m = 1,4



                         if ((trim(rdline)).EQ.string(m))  then

                                read(12,'(A)',iostat=filestat), arguments(m,1)

                                ! print *, arguments(m,1)
                                exit
                        end if



                end do

         end if


 !!! MISCELLANEOUS STUFF

 end do readline
 end if


END SUBROUTINE  Basis_parsing_subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!!!!!!!

SUBROUTINE  Region_parsing_subroutine(string,arguments,rdline)

!!!! initializing values
 Implicit None

 character(len=100)                         :: arguments(:,:)
 character(len=100)                         :: rdline
 character(len=100)                         :: string(:)
 integer                                    :: m

 if (trim(rdline)=="START_REGION") then

 readline:do

        read(12,'(A)',iostat=filestat), rdline
        ! terminate upon encoutering END_start_string
        if (trim(rdline)=="END_REGION" ) then

                exit

        end if

! to remove spaces

        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")

 !!!! TO IGNORE BLANK SPACES AND COMMENTS

        if ( (trim(rdline).NE."") .AND.  (trim(rdline(:1)).NE."!")) then
!!! search for key words and store values


                  do m = 1,1


                         if ((trim(rdline)).EQ.string(m))  then

                                read(12,'(A)',iostat=filestat), arguments(m,1)

                                ! print *, arguments(m,1)
                                exit
                        end if



                end do

         end if


 !!! MISCELLANEOUS STUFF

 end do readline
 end if


END SUBROUTINE  Region_parsing_subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  CoordinateSystem_parsing_subroutine(string,arguments,rdline)

!!!! initializing values
 Implicit None

 character(len=100)                         :: arguments(:,:)
 character(len=100)                         :: rdline
 character(len=100)                         :: string(:)
 integer                                    :: m

! print *, trim(rdline)
! print *, trim(rdline)=="START_COORDINATE_SYSTEM"
 if (trim(rdline)=="START_COORDINATE_SYSTEM") then

 readline:do

        read(12,'(A)',iostat=filestat), rdline
        ! terminate upon encoutering END_start_string
        if (trim(rdline)=="END_COORDINATE_SYSTEM" ) then

                exit

        end if

! to remove spaces

        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")

 !!!! TO IGNORE BLANK SPACES AND COMMENTS

        if ( (trim(rdline).NE."") .AND.  (trim(rdline(:1)).NE."!")) then
!!! search for key words and store values


                  do m = 1,2


                         if ((trim(rdline)).EQ.string(m))  then

                                read(12,'(A)',iostat=filestat), arguments(m,1)

                                ! print *, arguments(m,1)
                                exit
                        end if



                end do

         end if


 !!! MISCELLANEOUS STUFF

 end do readline
 end if


END SUBROUTINE  CoordinateSystem_parsing_subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  ControlLoop_parsing_subroutine(string,arguments,rdline)

!!!! initializing values
 Implicit None

 character(len=100)                         :: arguments(:,:)
 character(len=100)                         :: rdline
 character(len=100)                         :: string(:)
 integer                                    :: m


 if (trim(rdline)=="START_CONTROL_LOOP") then

 readline:do

        read(12,'(A)',iostat=filestat), rdline
        ! terminate upon encoutering END_start_string
        if (trim(rdline)=="END_CONTROL_LOOP" ) then

                exit

        end if

! to remove spaces

        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")

 !!!! TO IGNORE BLANK SPACES AND COMMENTS

        if ( (trim(rdline).NE."") .AND.  (trim(rdline(:1)).NE."!")) then
!!! search for key words and store values


                  do m = 1,3


                         if ((trim(rdline)).EQ.string(m))  then

                                read(12,'(A)',iostat=filestat), arguments(m,1)

                                ! print *, arguments(m,1)
                                exit
                        end if



                end do

         end if


 !!! MISCELLANEOUS STUFF

 end do readline
 end if


END SUBROUTINE  ControlLoop_parsing_subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  DependentField_parsing_subroutine(string,arguments,rdline)

!!!! initializing values
 Implicit None

 character(len=100)                         :: arguments(:,:)
 character(len=100)                         :: rdline
 character(len=100)                         :: string(:)
 integer                                    :: m

 if (trim(rdline)=="START_DEPENDENT_FIELD") then

 readline:do

        read(12,'(A)',iostat=filestat), rdline
        ! terminate upon encoutering END_start_string
        if (trim(rdline)=="END_DEPENDENT_FIELD" ) then

                exit

        end if

! to remove spaces

        rdline =  strip_space(rdline)
        call remove_character(rdline,"!")

 !!!! TO IGNORE BLANK SPACES AND COMMENTS

        if ( (trim(rdline).NE."") .AND.  (trim(rdline(:1)).NE."!")) then
!!! search for key words and store values


                  do m = 1,2



                         if ((trim(rdline)).EQ.string(m))  then

                                read(12,'(A)',iostat=filestat), arguments(m,1)

                                ! print *, arguments(m,1)
                                exit
                        end if



                end do

         end if


 !!! MISCELLANEOUS STUFF

 end do readline
 end if


END SUBROUTINE  DependentField_parsing_subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE  searching(rdline,search_string,counter)

!!!! initializing values
 Implicit None


 character(len=100)                         :: rdline
 character(len=*)                           :: search_string
 integer                                    :: counter

rdline =  strip_space(rdline)

!search_string = strip_space(search_string)
!! print *, rdline
! search_string
!! print *, trim(rdline)==search_string
if (trim(rdline)==search_string) then
    counter = counter+1
!    ! print *, "yes i did it , did it again , i palyed with your heart"
end if

END SUBROUTINE  searching



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








FUNCTION strip_space( input_string ) RESULT ( output_string )

    ! -- Arguments
    CHARACTER( * ), INTENT( IN )  :: input_string
    INTEGER                       :: n

    ! -- Function result
    CHARACTER( LEN( input_string ) ) :: output_string

    ! -- Local parameters
    INTEGER,        PARAMETER :: IACHAR_SPACE = 32, &
                                 IACHAR_TAB   = 9

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

      ! -- If the character is NOT a space ' ' or a tab '->|'
      ! -- copy it to the output string.
      IF ( iachar_character /= IACHAR_SPACE .AND. &
           iachar_character /= IACHAR_TAB         ) THEN
        n = n + 1
        output_string( n:n ) = input_string( i:i )
      END IF

    END DO


  END FUNCTION strip_space


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.REMOVE A CHARACTER FROM A STRING!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine remove_character(string,character)
  implicit none
  integer :: ind
  character(len=100)  :: string
  character           :: character
  ind = SCAN(string, character)
  if (ind .NE.0) then
  string = string(:ind-1)
  end if
end subroutine remove_character




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! STORE AND SPLITS ARGUMENTS !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine store_argument(matrix)
  character(len=100)                 :: string                     ! initializing variables
  character(len=100),   dimension(:) :: matrix
  integer :: g,k,filestat
  g=0
  !!!! TERMINATE UPON ENDIF

  g=g+1
  read(string,*) matrix(g)
  k=1
  ! print *, matrix(g)

end subroutine store_argument



 subroutine split_string(str,word)

  character(len = 100)                   :: str
  CHARACTER(len = 10),dimension(*)       :: word
  INTEGER                                :: pos1 = 1, pos2, n = 0, i

   DO
    pos2 = INDEX(str(pos1:), ",")
    IF (pos2 == 0) THEN
       n = n + 1
       word(n) = str(pos1:)
       EXIT
    END IF
    n = n + 1
    word(n) = str(pos1:pos1+pos2-2)
    pos1 = pos2+pos1

   END DO

 end subroutine split_string


REAL(CMISSRP) function str2real(str_2)
    implicit none
    ! Arguments
    character(len=*)      :: str_2
    read(str_2,*,iostat=stat)  str2real
end function str2real

integer function str2int(str)
    implicit none
    ! Arguments
    character(len=*)      :: str
    read(str,*,iostat=stat)  str2int
end function str2int


end module parsing
