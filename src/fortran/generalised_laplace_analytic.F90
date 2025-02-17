PROGRAM GeneralisedLaplaceAnalytic

  USE OpenCMISS
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Program parameters
  REAL(OC_RP), PARAMETER :: WIDTH=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
 
  REAL(OC_RP), PARAMETER :: SIGMA_11=3.0_OC_RP
  REAL(OC_RP), PARAMETER :: SIGMA_22=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: SIGMA_12=0.0_OC_RP  
  
  REAL(OC_RP), PARAMETER :: PI=3.141592653589793_OC_RP

  INTEGER(OC_Intg), PARAMETER :: MAX_NUMBER_OF_REFINEMENTS=5
  INTEGER(OC_Intg), PARAMETER :: MAX_NUMBER_OF_INTERPOLATIONS=9
  
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=7
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=8
  INTEGER(OC_Intg), PARAMETER :: FIBRE_FIELD_USER_NUMBER=9
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program types

  !Program variables
  INTEGER(OC_Intg) :: argumentLength,numberOfArguments,status
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex
  INTEGER(OC_Intg) :: baseNumberOfGlobalXElements,baseNumberOfGlobalYElements,fibreAngleDegrees,interpolationType, &
    & numberOfGaussXi,numberOfGlobalXElements,numberOfGlobalYElements
  INTEGER(OC_Intg) :: err
  REAL(OC_RP) :: fibreAngleRadians
  INTEGER(OC_Intg) :: convergenceDataIdx,numberOfDOFs,numberOfInterpolations,numberOfRefinements,refinementIdx
  REAL(OC_RP) :: integralValues(2),localValues(8),localGhostValues(8),ghostIntegralValues(2), globalValues(8)
  REAL(OC_RP) :: convergenceData(4,MAX_NUMBER_OF_REFINEMENTS,MAX_NUMBER_OF_INTERPOLATIONS), &
    & logConvergenceData(4,MAX_NUMBER_OF_REFINEMENTS,MAX_NUMBER_OF_INTERPOLATIONS), &
    & meanLogConvergenceData(4,MAX_NUMBER_OF_INTERPOLATIONS),slopeConvergenceData(4,MAX_NUMBER_OF_INTERPOLATIONS), &
    & sumLogConvergenceData(6,MAX_NUMBER_OF_INTERPOLATIONS)
  CHARACTER(LEN=255) :: commandArgument,filename

  !OpenCMISS variables
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: analyticField,dependentField,geometricField,equationsSetField,fibreField,materialsField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_NodesType) :: nodes
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 4) THEN
    !If we have enough arguments then use the first three for setting up the problem.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) baseNumberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of base X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) baseNumberOfGlobalYElements
    IF(numberOfGlobalYElements<=0) CALL HandleError("Invalid number of base Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) fibreAngleDegrees
    IF(fibreAngleDegrees<=-180_OC_Intg.Or.fibreAngleDegrees>180_OC_Intg) CALL HandleError("Invalid fibre angle in degrees.")
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) numberOfRefinements
    IF(numberOfRefinements<=0) CALL HandleError("Invalid number of refinements.")
    IF(numberOfRefinements>MAX_NUMBER_OF_REFINEMENTS) CALL HandleError("Increase the maximum number of refinements.")
  ELSE
    !If there are not enough arguments default the problem specification
    baseNumberOfGlobalXElements=4
    baseNumberOfGlobalYElements=2
    fibreAngleDegrees=60_OC_Intg
    numberOfRefinements=4
  ENDIF
  fibreAngleRadians=REAL(fibreAngleDegrees,OC_RP)*2.0_OC_RP*PI/360.0_OC_RP

  numberOfInterpolations=3

  convergenceData=0.0_OC_RP
  logConvergenceData=0.0_OC_RP
  meanLogConvergenceData=0.0_OC_RP
  slopeConvergenceData=0.0_OC_RP
  sumLogConvergenceData=0.0_OC_RP

  !Initialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  
  DO interpolationType=1,numberOfInterpolations
    
    DO refinementIdx=1,numberOfRefinements

      numberOfGlobalXElements=refinementIdx*baseNumberOfGlobalXElements
      numberOfGlobalYElements=refinementIdx*baseNumberOfGlobalYElements

      convergenceData(1,refinementIdx,interpolationType)=SQRT((WIDTH*HEIGHT)/(REAL(numberOfGlobalXElements,OC_RP)* &
        & REAL(numberOfGlobalYElements,OC_RP)))

      !Create the context
      CALL OC_Context_Initialise(context,err)
      CALL OC_Context_Create(refinementIdx,context,err)
      CALL OC_Region_Initialise(worldRegion,err)
      CALL OC_Context_WorldRegionGet(context,worldRegion,err)
      CALL OC_Context_RandomSeedsSet(context,9999,err)
      !CALL OC_DiagnosticsSetOn(OC_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["Laplace_FiniteElementCalculate"],err)

      WRITE(filename,'(A,"_",I0,"x",I0,"_",I0,"_",I0)') "GeneralisedLaplace",numberOfGlobalXElements,numberOfGlobalYElements, &
        & interpolationType,fibreAngleDegrees

      CALL OC_OutputSetOn(filename,err)

      !Get the computational nodes information
      CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
      CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)

      CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
      CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
      CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
      CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

      !-----------------------------------------------------------------------------------------------------------
      ! COORDINATE SYSTEM
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a new RC coordinate system
      CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
      CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
      !Set the coordinate system to be 2D
      CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
      !Finish the creation of the coordinate system
      CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! REGION
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the region
      CALL OC_Region_Initialise(region,err)
      CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
      !Set the regions coordinate system to the 2D RC coordinate system that we have created
      CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
      !Set the region label
      CALL OC_Region_LabelSet(region,"GeneralisedLaplace",err)
      !Finish the creation of the region
      CALL OC_Region_CreateFinish(region,err)

      !-----------------------------------------------------------------------------------------------------------
      ! BASIS
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a basis (default is trilinear lagrange)
      CALL OC_Basis_Initialise(basis,err)
      CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
      SELECT CASE(interpolationType)
      CASE(1,2,3,4)
        CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
      CASE(7,8,9)
        CALL OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE,err)
      CASE DEFAULT
        CALL HandleError("Invalid interpolation type.")
      END SELECT
      SELECT CASE(interpolationType)
      CASE(1)
        numberOfGaussXi=2
      CASE(2)
        numberOfGaussXi=3
      CASE(3,4)
        numberOfGaussXi=4
      CASE DEFAULT
        numberOfGaussXi=0 !Don't set number of Gauss points for tri/tet
      END SELECT
      !Set the basis to be a bi-interpolation basis
      CALL OC_Basis_NumberOfXiSet(basis,2,err)
      CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
      IF(numberOfGaussXi>0) THEN
        CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
      ENDIF
      !Finish the creation of the basis
      CALL OC_Basis_CreateFinish(basis,err)

      !-----------------------------------------------------------------------------------------------------------
      ! MESH
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a generated mesh in the region
      CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
      CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
      !Set up a regular x*y*z mesh
      CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
      !Set the default basis
      CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)
      !Define the mesh on the region
      CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
      CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
      !Finish the creation of a generated mesh in the region
      CALL OC_Mesh_Initialise(mesh,err)
      CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

      !Create a decomposition
      CALL OC_Decomposition_Initialise(decomposition,err)
      CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
      !Finish the decomposition
      CALL OC_Decomposition_CreateFinish(decomposition,err)

      !Decompose
      CALL OC_Decomposer_Initialise(decomposer,err)
      CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
      !Add in the decomposition
      CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
      !Finish the decomposer
      CALL OC_Decomposer_CreateFinish(decomposer,err)

      !Destory the mesh now that we have decomposed it
      !CALL OC_Mesh_Destroy(mesh,err)

      !-----------------------------------------------------------------------------------------------------------
      ! GEOMETRIC FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Start to create a default (geometric) field on the region
      CALL OC_Field_Initialise(geometricField,err)
      CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
      !Set the decomposition to use
      CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
      !Set the variable label
      CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
      !Set the domain to be used by the field components.
      CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
      CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,1,err)
      !Finish creating the field
      CALL OC_Field_CreateFinish(geometricField,err)

      !Update the geometric field parameters
      CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

      !-----------------------------------------------------------------------------------------------------------
      ! FIBRE FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Start to create a fibre field on the region
      CALL OC_Field_Initialise(fibreField,err)
      CALL OC_Field_CreateStart(FIBRE_FIELD_USER_NUMBER,region,fibreField,err)
      !Set the field type
      CALL OC_Field_TypeSet(fibreField,OC_FIELD_FIBRE_TYPE,err)
      !Set the decomposition to use
      CALL OC_Field_DecompositionSet(fibreField,decomposition,err)
      !Set the geometric field
      CALL OC_Field_GeometricFieldSet(fibreField,geometricField,err)
      !Set the variable label
      CALL OC_Field_VariableLabelSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",err)
      !Set the domain to be used by the field components.
      CALL OC_Field_ComponentMeshComponentSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
      !Set the interpolation type to be constant
      CALL OC_Field_ComponentInterpolationSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_CONSTANT_INTERPOLATION,err)
      CALL OC_Field_ComponentInterpolationSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,2,OC_FIELD_CONSTANT_INTERPOLATION,err)
      !Finish creating the field
      CALL OC_Field_CreateFinish(fibreField,err)

      !Initialise the field parameters to the fibre angle (in radians)
      CALL OC_Field_ComponentValuesInitialise(fibreField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
        & fibreAngleRadians,err)

      !-----------------------------------------------------------------------------------------------------------
      ! EQUATIONS SETS
      !-----------------------------------------------------------------------------------------------------------

      !Create the Generalised Laplace equations set
      CALL OC_EquationsSet_Initialise(equationsSet,err)
      CALL OC_Field_Initialise(equationsSetField,err)
      CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
        & OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,OC_EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE], &
        & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
      !Finish creating the equations set
      CALL OC_EquationsSet_CreateFinish(equationsSet,err)

      !-----------------------------------------------------------------------------------------------------------
      ! DEPENDENT FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set dependent field variables
      CALL OC_Field_Initialise(dependentField,err)
      CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
      !Finish the equations set dependent field variables
      CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

      !Initialise the field with an initial guess
      CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
        & 0.0_OC_RP, err)

      !Get the number of DOFs
      CALL OC_Field_NumberOfGlobalDOFsGet(dependentField,OC_FIELD_U_VARIABLE_TYPE,numberOfDOFs,err)
      convergenceData(2,refinementIdx,interpolationType)=REAL(numberOfDOFs,OC_RP)

      !-----------------------------------------------------------------------------------------------------------
      ! MATERIALS FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set materials field variables
      CALL OC_Field_Initialise(materialsField,err)
      CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
      !Set the interpolation type to be constant
      CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,1,OC_FIELD_CONSTANT_INTERPOLATION,err)
      CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,2,OC_FIELD_CONSTANT_INTERPOLATION,err)
      CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,3,OC_FIELD_CONSTANT_INTERPOLATION,err)
      !Finish the equations set materials field variables
      CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

      !Initialise the materials conductivity tensor
      CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,SIGMA_11,err)
      CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,SIGMA_22,err)
      CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,SIGMA_12,err)

      !-----------------------------------------------------------------------------------------------------------
      ! ANALYTIC FIELD
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set analytic field variables
      CALL OC_Field_Initialise(analyticField,err)
      CALL OC_EquationsSet_AnalyticCreateStart(equationsSet,OC_EQUATIONS_SET_GENERALISED_LAPLACE_EQUATION_TWO_DIM_1, &
        & ANALYTIC_FIELD_USER_NUMBER,analyticField,Err)
      !Finish the equations set analytic field variables
      CALL OC_EquationsSet_AnalyticCreateFinish(equationsSet,err)

      !Don't set the analytic field values, they will default to the materials and fibre field values

      !-----------------------------------------------------------------------------------------------------------
      ! EQUATIONS
      !-----------------------------------------------------------------------------------------------------------

      !Create the equations set equations
      CALL OC_Equations_Initialise(equations,err)
      CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
      !Set the equations matrices sparsity type
      CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
      !CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_FULL_MATRICES,err)
      !Set the equations set output
      !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
      !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_TIMING_OUTPUT,err)
      !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
      CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
      !Finish the equations set equations
      CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

      !-----------------------------------------------------------------------------------------------------------
      ! PROBLEM
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of a problem.
      CALL OC_Problem_Initialise(problem,err)
      CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
        & OC_PROBLEM_LAPLACE_EQUATION_TYPE,OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE],problem,err)
      !Finish the creation of a problem.
      CALL OC_Problem_CreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! CONTROL LOOP
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the problem control loop
      CALL OC_Problem_ControlLoopCreateStart(problem,err)
      !Finish creating the problem control loop
      CALL OC_Problem_ControlLoopCreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! SOLVER
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the problem solvers
      CALL OC_Solver_Initialise(solver,err)
      CALL OC_Problem_SolversCreateStart(problem,err)
      CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
      CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_NO_OUTPUT,err)
      !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
      !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT,err)
      !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT,err)
      !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)

      CALL OC_Solver_LinearTypeSet(solver,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
      CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(solver,1.0E-12_OC_RP,err)
      CALL OC_Solver_LinearIterativeRelativeToleranceSet(solver,1.0E-12_OC_RP,err)

      !CALL OC_Solver_LinearTypeSet(solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
      !CALL OC_Solver_LinearTypeSet(solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
      !CALL OC_Solver_LibraryTypeSet(solver,OC_SOLVER_MUMPS_LIBRARY,err)
      !CALL OC_Solver_LibraryTypeSet(solver,OC_SOLVER_LAPACK_LIBRARY,err)
      !CALL OC_Solver_LibraryTypeSet(solver,OC_SOLVER_SUPERLU_LIBRARY,err)
      !CALL OC_Solver_LibraryTypeSet(solver,OC_SOLVER_PASTIX_LIBRARY,err)
      !Finish the creation of the problem solver
      CALL OC_Problem_SolversCreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! SOLVER EQUATIONS
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the problem solver equations
      CALL OC_Solver_Initialise(solver,err)
      CALL OC_SolverEquations_Initialise(solverEquations,err)
      CALL OC_Problem_SolverEquationsCreateStart(problem,err)
      !Get the solve equations
      CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
      CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
      !Set the solver equations sparsity
      CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
      !CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES,err)
      !Add in the equations set
      CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
      !Finish the creation of the problem solver equations
      CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! BOUNDARY CONDITIONS
      !-----------------------------------------------------------------------------------------------------------

      !Start the creation of the equations set boundary conditions
      CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
      CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
      !Set the analytic boundary conditions
      CALL OC_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
      !Finish the creation of the equations set boundary conditions
      CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

      !-----------------------------------------------------------------------------------------------------------
      ! SOLVE
      !-----------------------------------------------------------------------------------------------------------

      !Solve the problem
      CALL OC_Problem_Solve(problem,err)

      !-----------------------------------------------------------------------------------------------------------
      ! OUTPUT
      !-----------------------------------------------------------------------------------------------------------

      !Perform Analytic analysis
      CALL OC_AnalyticAnalysis_Output(dependentField,filename,err)
      CALL OC_AnalyticAnalysis_RMSErrorGetNode(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,OC_ANALYTIC_ABSOLUTE_ERROR_TYPE, &
        & localValues,localGhostValues,globalValues,err)
      convergenceData(3,refinementIdx,interpolationType)=globalValues(1)
      CALL OC_AnalyticAnalysis_IntegralNumericalValueGet(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,integralValues,&
        & ghostIntegralValues,err)
      convergenceData(3,refinementIdx,interpolationType)=integralValues(1)

      !Export results
      CALL OC_Fields_Initialise(fields,err)
      CALL OC_Fields_Create(region,fields,err)
      CALL OC_Fields_NodesExport(fields,filename,"FORTRAN",err)
      CALL OC_Fields_ElementsExport(fields,filename,"FORTRAN",err)
      CALL OC_Fields_Finalise(fields,err)

      !Turn ouput off
      CALL OC_OutputSetOff(err)
      
      !Destroy the context
      CALL OC_Context_Destroy(context,err)
      
      DO convergenceDataIdx=1,4
        logConvergenceData(convergenceDataIdx,refinementIdx,interpolationType)= &
          & LOG10(convergenceData(convergenceDataIdx,refinementIdx,interpolationType))
        meanLogConvergenceData(convergenceDataIdx,refinementIdx)=meanLogConvergenceData(convergenceDataIdx,refinementIdx)+ &
          & logConvergenceData(convergenceDataIdx,refinementIdx,interpolationType)
      ENDDO !convergenceDataIdx
      
    ENDDO !refinementIdx

    !Calculate the slope of the convergence line
    meanLogConvergenceData(:,interpolationType)=meanLogConvergenceData(:,interpolationType)/REAL(numberOfRefinements,OC_RP)
    DO refinementIdx=1,numberOfRefinements
      sumLogConvergenceData(1,interpolationType)=sumLogConvergenceData(1,interpolationType)+ &
        & (logConvergenceData(1,refinementIdx,interpolationType)-meanLogConvergenceData(1,interpolationType))**2
      sumLogConvergenceData(2,interpolationType)=sumLogConvergenceData(2,interpolationType)+ &
        & (logConvergenceData(2,refinementIdx,interpolationType)-meanLogConvergenceData(2,interpolationType))**2
      sumLogConvergenceData(3,interpolationType)=sumLogConvergenceData(3,interpolationType)+ &
        & (logConvergenceData(1,refinementIdx,interpolationType)-meanLogConvergenceData(1,interpolationType))* &
        & (logConvergenceData(3,refinementIdx,interpolationType)-meanLogConvergenceData(3,interpolationType))
      sumLogConvergenceData(4,interpolationType)=sumLogConvergenceData(4,interpolationType)+ &
        & (logConvergenceData(1,refinementIdx,interpolationType)-meanLogConvergenceData(1,interpolationType))* &
        & (logConvergenceData(4,refinementIdx,interpolationType)-meanLogConvergenceData(4,interpolationType))
      sumLogConvergenceData(5,interpolationType)=sumLogConvergenceData(5,interpolationType)+ &
        & (logConvergenceData(2,refinementIdx,interpolationType)-meanLogConvergenceData(2,interpolationType))* &
        & (logConvergenceData(3,refinementIdx,interpolationType)-meanLogConvergenceData(3,interpolationType))
      sumLogConvergenceData(6,interpolationType)=sumLogConvergenceData(6,interpolationType)+ &
        & (logConvergenceData(2,refinementIdx,interpolationType)-meanLogConvergenceData(2,interpolationType))* &
        & (logConvergenceData(4,refinementIdx,interpolationType)-meanLogConvergenceData(4,interpolationType))
    ENDDO ! refinementIdx
    slopeConvergenceData(1,interpolationType)=sumLogConvergenceData(3,interpolationType)/sumLogConvergenceData(1,interpolationType)
    slopeConvergenceData(2,interpolationType)=sumLogConvergenceData(4,interpolationType)/sumLogConvergenceData(1,interpolationType)
    slopeConvergenceData(3,interpolationType)=sumLogConvergenceData(5,interpolationType)/sumLogConvergenceData(2,interpolationType)
    slopeConvergenceData(4,interpolationType)=sumLogConvergenceData(6,interpolationType)/sumLogConvergenceData(2,interpolationType)

    WRITE(*,'(X)')
    WRITE(*,'("Interpolation type : ",I0)') interpolationType
    WRITE(*,'("  RMS error vs h slope            = ",E12.4)') slopeConvergenceData(1,interpolationType)
    WRITE(*,'("  Integral^2 error vs h slope     = ",E12.4)') slopeConvergenceData(2,interpolationType)
    WRITE(*,'("  RMS error vs #DOFs slope        = ",E12.4)') slopeConvergenceData(3,interpolationType)
    WRITE(*,'("  Integral^2 error vs #DOFs slope = ",E12.4)') slopeConvergenceData(4,interpolationType)

  ENDDO !interpolationType
  
  !Finialise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."
  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM GeneralisedLaplaceAnalytic
