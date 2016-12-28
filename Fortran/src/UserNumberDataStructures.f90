
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

!!!!!!!!!!!!!!!!!! Assigning tags. NOTE:  Each tag should be distinct !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
j = 1  

do i = 1, num_of_CoordinateSystem

  	CoordinateSystemUserNumber(num_of_CoordinateSystem)  = j 
        j = j + 1 
end do 

do i = 1, num_of_Region

  	RegionUserNumber(num_of_Region)  = j 
        j = j + 1 
end do 


do i = 1, num_of_Basis

  	BasisUserNumber(num_of_Basis)  = j 
        j = j + 1 
end do 


do i = 1, num_of_PressureBasis

  	PressureBasisUserNumber(num_of_PressureBasis)  = j 
        j = j + 1 
end do 

do i = 1,num_of_GeneratedMesh

  	GeneratedMeshUserNumber(num_of_GeneratedMesh)  = j 
        j = j + 1 
end do 

do i = 1, num_of_Mesh

 	MeshUserNumber(num_of_Mesh) = j 
        j = j + 1 
end do 

do i = 1, num_of_Decomposition

 	DecompositionUserNumber(num_of_Decomposition) = j 
        j = j + 1 
end do 

do i = 1, num_of_GeometricField

 	FieldGeometryUserNumber(num_of_GeometricField) = j 
        j = j + 1 
end do

do i = 1, num_of_FiberField

 	FieldFibreUserNumber(num_of_FiberField) = j 
        j = j + 1 
end do

do i = 1, num_of_MaterialField

 	FieldMaterialUserNumber(num_of_MaterialField) = j 
        j = j + 1 
end do

do i = 1, num_of_DependentField

 	FieldDependentUserNumber(num_of_DependentField) = j 
        j = j + 1 
end do

do i = 1, num_of_EquationsSet

 	EquationSetUserNumber(num_of_EquationsSet) = j 
        j = j + 1 
end do

do i = 1, num_of_EquationSetField

 	EquationsSetFieldUserNumber(num_of_EquationSetField) = j 
        j = j + 1 
end do


do i = 1, num_of_Problem

 	ProblemUserNumber(num_of_Problem) = j 
        j = j + 1 
end do
