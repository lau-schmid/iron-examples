

 !!!!!!!!!!!!!!!!!!!!________THE FOLLOWING HEADER FILE ALLOCATES THE SIZE OF "USER NUMBER" DATA STRUCTURES BASED ON THE VARIABLES INITIALIED IN _____!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!! ....... HEADER FILE ALLOCATING DERIVED DATA STRUCTURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________________________________________________________________!!!!!!!!!!!!!!!!

  allocate (CoordinateSystemUserNumber(NumberOfCoordinateSystem))
  allocate (RegionUserNumber(NumberOfRegion))
  allocate (BasisUserNumber(NumberOfBasis))
  allocate (PressureBasisUserNumber(NumberOfPressureBasis))
  allocate (GeneratedMeshUserNumber(NumberOfGeneratedMesh))
  allocate (MeshUserNumber(NumberOfMesh))
  allocate (DecompositionUserNumber(NumberOfDecomposition))
  allocate (FieldGeometryUserNumber(NumberOfGeometricField))
  allocate (FieldFibreUserNumber(NumberOfFiberField))
  allocate (FieldMaterialUserNumber(NumberOfMaterialField))
  allocate (FieldDependentUserNumber(NumberOfDependentField))
  allocate (EquationSetUserNumber(NumberOfEquationsSet))
  allocate (EquationsSetFieldUserNumber(NumberOfEquationSetField ))
  allocate (ProblemUserNumber(NumberOfProblem))
  allocate (FieldUserNumber(NumberOfProblem))

!!!!!!!!!!!!!!!!!! Assigning tags. NOTE:  Each tag should be distinct !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 UserNumberLabel = 1

  do i = 1, NumberOfCoordinateSystem

        CoordinateSystemUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
  end do

 do i = 1, NumberOfRegion

        RegionUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 end do


 do i = 1, NumberOfBasis+1

        BasisUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 end do


 do i = 1, NumberOfPressureBasis

        PressureBasisUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1,NumberOfGeneratedMesh

        GeneratedMeshUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfMesh

       MeshUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfDecomposition

       DecompositionUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfGeometricField

       FieldGeometryUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfFiberField

       FieldFibreUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfMaterialField

       FieldMaterialUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfDependentField

       FieldDependentUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfEquationsSet

       EquationSetUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfEquationsSet

       EquationsSetFieldUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfProblem

       ProblemUserNumber(i) = UserNumberLabel
       UserNumberLabel = UserNumberLabel + 1
 end do

 do i = 1, NumberOfField

        FieldUserNumber(i)  = UserNumberLabel
        UserNumberLabel = UserNumberLabel + 1
 end do

