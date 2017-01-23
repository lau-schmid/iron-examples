
  allocate (CoordinateSystemUserNumber(num_of_CoordinateSystem))  
  allocate (RegionUserNumber(num_of_Region))
  allocate (BasisUserNumber(num_of_Basis))
  allocate (PressureBasisUserNumber(num_of_PressureBasis))
  allocate (GeneratedMeshUserNumber(num_of_GeneratedMesh))
  allocate (MeshUserNumber(num_of_Mesh))
  allocate (DecompositionUserNumber(num_of_Decomposition))
  allocate (FieldGeometryUserNumber(num_of_GeometricField))
  allocate (FieldFibreUserNumber(num_of_FiberField))
  allocate (FieldMaterialUserNumber(num_of_MaterialField))
  allocate (FieldDependentUserNumber(num_of_DependentField))
  allocate (EquationSetUserNumber(num_of_EquationsSet))
  allocate (EquationsSetFieldUserNumber(num_of_EquationSetField ))
  allocate (ProblemUserNumber(num_of_Problem))
  allocate (FieldUserNumber(num_of_Field))
!!!!!!!!!!!!!!!!!! Assigning tags. NOTE:  Each tag should be distinct !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
j = 1  

do i = 1, num_of_CoordinateSystem

  	CoordinateSystemUserNumber(i)  = j 
        j = j + 1 
end do 

do i = 1, num_of_Region

  	RegionUserNumber(i)  = j 
        j = j + 1 
end do 


do i = 1, num_of_Basis+1

  	BasisUserNumber(i)  = j 
        j = j + 1 
end do 


do i = 1, num_of_PressureBasis

  	PressureBasisUserNumber(i)  = j 
        j = j + 1 
end do 

do i = 1,num_of_GeneratedMesh

  	GeneratedMeshUserNumber(i)  = j 
        j = j + 1 
end do 

do i = 1, num_of_Mesh

 	MeshUserNumber(i) = j 
        j = j + 1 
end do 

do i = 1, num_of_Decomposition

 	DecompositionUserNumber(i) = j 
        j = j + 1 
end do 

do i = 1, num_of_GeometricField

 	FieldGeometryUserNumber(i) = j 
        j = j + 1 
end do

do i = 1, num_of_FiberField

 	FieldFibreUserNumber(i) = j 
        j = j + 1 
end do

do i = 1, num_of_MaterialField

 	FieldMaterialUserNumber(i) = j 
        j = j + 1 
end do

do i = 1, num_of_DependentField

 	FieldDependentUserNumber(i) = j 
        j = j + 1 
end do

do i = 1, num_of_EquationsSet

 	EquationSetUserNumber(i) = j 
        j = j + 1 
end do

do i = 1, 1

 	EquationsSetFieldUserNumber(i) = j 

        j = j + 1 
end do
        

do i = 1, num_of_Problem

 	ProblemUserNumber(i) = j 
        j = j + 1 
end do

do i = 1, num_of_Field

  	FieldUserNumber(i)  = j 
        j = j + 1 
end do 
