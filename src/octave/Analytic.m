# Constant declarations

global DIAGNOSTICS_OFF = 0;
global DIAGNOSTICS_ON = 1;

global LINEAR_LAGRANGE_INTERPOLATION = 1;
global QUADRATIC_LAGRANGE_INTERPOLATION = 2;
global CUBIC_LAGRANGE_INTERPOLATION = 3;

global OFF_BOUNDARY_NODE = 0;
global ON_BOUNDARY_NODE = 1;

global DOF_FIXED = 0;
global DOF_FREE = 1;

global GENERALISED_LAPLACE_ANALYTIC_TYPE = 1;
global STANDARD_LAPLACE_ANALYTIC_TYPE = 2;

# Example parameters and options

global diagnostics = DIAGNOSTICS_OFF;
#global diagnostics = DIAGNOSTICS_ON;

#global analyticType = STANDARD_LAPLACE_ANALYTIC_TYPE;
global analyticType = GENERALISED_LAPLACE_ANALYTIC_TYPE;

#global L = 1.0;
global L = 2.0;
global H = 1.0;

global sigmat = 3.0;
global sigman = 1.0;
#global thetadeg = 30.0;
global thetadeg = 60.0;
#global sigmat = 1.0;
#global sigman = 1.0;
#global thetadeg = 0.0;
global theta = thetadeg*pi()/180.0;

#global baseNumXElements = 3;
global baseNumXElements = 4;
global baseNumYElements = 2;
#global baseNumYElements = 3;
global refinementLevel = 2;

interpolationType = LINEAR_LAGRANGE_INTERPOLATION;

# Should not need to change below here

global sigma11 = sigmat*cos(theta)*cos(theta)+sigman*sin(theta)*sin(theta);
global sigma12 = 1.0*(sigmat-sigman)*sin(theta)*cos(theta);
global sigma21 = 1.0*(sigmat-sigman)*sin(theta)*cos(theta);
global sigma22 = sigmat*sin(theta)*sin(theta)+sigman*cos(theta)*cos(theta);
global sigma = [ sigma11, sigma12 ; sigma21, sigma22 ];
global lambda1 = sigma12/sigma22;
global lambda2 = sqrt(sigma11*sigma22-sigma12*sigma21)/sigma22;

global numXElements = baseNumXElements * refinementLevel
global numYElements = baseNumYElements * refinementLevel
global numElements = numXElements * numYElements
global numNodesPerXi = (interpolationType + 1)
global numElementDOFS = numNodesPerXi*numNodesPerXi
global numXNodes = (numNodesPerXi-1)*numXElements + 1
global numYNodes = (numNodesPerXi-1)*numYElements + 1
global numNodes = numXNodes*numYNodes
global numDOFS = numNodes
global numFixedDOFS = 2*numXNodes + 2*(numYNodes-2)
global numFreeDOFS = numDOFS - numFixedDOFS

function [ nodePositions, boundaryNodes ] = computeGeometry( interpolationType )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global ON_BOUNDARY_NODE;
  global OFF_BOUNDARY_NODE;
  global LINEAR_LAGRANGE_INTERPOLATION;
  global QUADRATIC_LAGRANGE_INTERPOLATION;
  global CUBIC_LAGRANGE_INTERPOLATION;
  global L;
  global H;
  global numXElements;
  global numYElements;
  global numNodes;

  nodeIdx = 0;
  for j = 1:numYElements+1
    for i = 1:numXElements+1
      nodeIdx = nodeIdx + 1;
      boundaryNodes(nodeIdx,1) = OFF_BOUNDARY_NODE;
      nodePositions(nodeIdx,1) = (i-1)*L/numXElements;
      nodePositions(nodeIdx,2) = (j-1)*H/numYElements;
      if( or( i==1, i==numXElements+1, j==1, j==numYElements+1 ) )
        boundaryNodes(nodeIdx,1) = ON_BOUNDARY_NODE;
      endif
    endfor
  endfor

  if( diagnostics == DIAGNOSTICS_ON )
    nodePositions
    boundaryNodes
  endif

endfunction

function [ elementDOFNumbers ] = computeElementDOFs( interpolationType, xElementIdx, yElementIdx )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global LINEAR_LAGRANGE_INTERPOLATION;
  global QUADRATIC_LAGRANGE_INTERPOLATION;
  global CUBIC_LAGRANGE_INTERPOLATION;
  global numXElements;
  global numYElements;

  switch( interpolationType )
    case LINEAR_LAGRANGE_INTERPOLATION
      elementDOFNumbers(1) = xElementIdx + (yElementIdx - 1)*(numXElements + 1);
      elementDOFNumbers(2) = 1 + xElementIdx + (yElementIdx - 1)*(numXElements + 1);
      elementDOFNumbers(3) = xElementIdx + yElementIdx*(numXElements + 1);
      elementDOFNumbers(4) = 1 + xElementIdx + yElementIdx*(numXElements + 1);
    case QUADRATIC_LAGRANGE_INTERPOLATION
      printf("ERROR: Quadratic Lagrange interpolation not implemented.\n");
      quit;
    case CUBIC_LAGRANGE_INTERPOLATION
      printf("ERROR: Cubic Lagrange interpolation not implemented.\n");
      quit;
    otherwise
      printf("ERROR: Invalid interpolation type.\n");
      quit;
  endswitch	  

  if( diagnostics == DIAGNOSTICS_ON )
    elementDOFNumbers
  endif
  
endfunction

function [ Ke, fe ] = computeElementStiffnessMatrix( interpolationType )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global LINEAR_LAGRANGE_INTERPOLATION;
  global QUADRATIC_LAGRANGE_INTERPOLATION;
  global CUBIC_LAGRANGE_INTERPOLATION;
  global L;
  global H;
  global numXElements;
  global numYElements;
  global sigma11;
  global sigma12;
  global sigma21;
  global sigma22;

  switch( interpolationType )
    case LINEAR_LAGRANGE_INTERPOLATION
      Ke(1,1) = (4.0*sigma11*H*H*numXElements*numXElements+3.0*(sigma12+sigma21)*L*H*numXElements*numYElements+4.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(1,2) = (-4.0*sigma11*H*H*numXElements*numXElements+3.0*(sigma12-sigma21)*L*H*numXElements*numYElements+2.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(1,3) = (2.0*sigma11*H*H*numXElements*numXElements-3.0*(sigma12-sigma21)*L*H*numXElements*numYElements-4.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(1,4) = (-2.0*sigma11*H*H*numXElements*numXElements-3.0*(sigma12+sigma21)*L*H*numXElements*numYElements-2.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(2,1) = Ke(1,2);
      Ke(2,2) = (4.0*sigma11*H*H*numXElements*numXElements-3.0*(sigma12+sigma21)*L*H*numXElements*numYElements+4.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(2,3) = (-2.0*sigma11*H*H*numXElements*numXElements+3.0*(sigma12+sigma21)*L*H*numXElements*numYElements-2.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(2,4) = (2.0*sigma11*H*H*numXElements*numXElements+3.0*(sigma12-sigma21)*L*H*numXElements*numYElements-4.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(3,1) = Ke(1,3);
      Ke(3,2) = Ke(2,3);
      Ke(3,3) = (4.0*sigma11*H*H*numXElements*numXElements-3.0*(sigma12+sigma21)*L*H*numXElements*numYElements+4.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(3,4) = (-4.0*sigma11*H*H*numXElements*numXElements-3.0*(sigma12-sigma21)*L*H*numXElements*numYElements+2.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      Ke(4,1) = Ke(1,4);
      Ke(4,2) = Ke(2,4);
      Ke(4,3) = Ke(3,4);
      Ke(4,4) = (4.0*sigma11*H*H*numXElements*numXElements+3.0*(sigma12+sigma21)*L*H*numXElements*numYElements+4.0*sigma22*L*L*numYElements*numYElements)/(12.0*L*H*numXElements*numYElements);
      fe(1,1) = 0.0;
      fe(2,1) = 0.0;
      fe(3,1) = 0.0;
      fe(4,1) = 0.0;
    case QUADRATIC_LAGRANGE_INTERPOLATION
      printf("ERROR: Quadratic Lagrange interpolation not implemented.\n");
      quit;
    case CUBIC_LAGRANGE_INTERPOLATION
      printf("ERROR: Cubic Lagrange interpolation not implemented.\n");
      quit;
    otherwise
      printf("ERROR: Invalid interpolation type.\n");
      quit;
  endswitch	  

  if( diagnostics == DIAGNOSTICS_ON )
    Ke
    fe
  endif
  
endfunction

function [ K, f ] = computeLinearMatrices( interpolationType )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global numXElements;
  global numYElements;
  global numElementDOFS;
  global numDOFS;

  K = zeros(numDOFS,numDOFS);
  f = zeros(numDOFS,1);

  Ke = zeros(numElementDOFS,numElementDOFS);
  fe = zeros(numElementDOFS,1);
  elementDOFNumbers = zeros(numElementDOFS,1);
  
  [ Ke, fe ] = computeElementStiffnessMatrix( interpolationType )

  elementIdx = 0;
  for yElementIdx = 1:numYElements
    for xElementIdx = 1:numXElements
      elementIdx=elementIdx + 1;
      [ elementDOFNumbers ] = computeElementDOFs( interpolationType, xElementIdx, yElementIdx );
      for rowDOFIdx = 1:numElementDOFS
	rowNumber = elementDOFNumbers( rowDOFIdx );
	for columnDOFIdx = 1:numElementDOFS
	  columnNumber = elementDOFNumbers( columnDOFIdx );
	  K(rowNumber,columnNumber) = K(rowNumber,columnNumber) + Ke(rowDOFIdx,columnDOFIdx);	
	endfor
	f(rowNumber,1) = f(rowNumber,1) + fe(rowDOFIdx,1);
      endfor
    endfor
  endfor
  
  if( diagnostics == DIAGNOSTICS_ON )
    K
    f
  endif
  
endfunction

function [ analyticU, gradAnalyticU, hessAnalyticU ] = analytic( analyticType, nodePositions )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global GENERALISED_LAPLACE_ANALYTIC_TYPE;
  global STANDARD_LAPLACE_ANALYTIC_TYPE;
  global lambda1;
  global lambda2;
  global numNodes;

  switch( analyticType )
    case GENERALISED_LAPLACE_ANALYTIC_TYPE
      for nodeIdx = 1:numNodes
	analyticU(nodeIdx,1) = 2.0*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*cos(lambda2*nodePositions(nodeIdx,2));
	gradAnalyticU(nodeIdx,1) = 2.0*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*cos(lambda2*nodePositions(nodeIdx,2));
	gradAnalyticU(nodeIdx,1) = -2.0*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*(lambda2*sin(lambda2*nodePositions(nodeIdx,2))+lambda1*cos(lambda2*nodePositions(nodeIdx,2)));
	hessAnalyticU(nodeIdx,1,1) = 2.0*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*cos(lambda2*nodePositions(nodeIdx,2));
	hessAnalyticU(nodeIdx,1,2) = -2.0*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*(lambda2*sin(lambda2*nodePositions(nodeIdx,2))+lambda1*cos(lambda2*nodePositions(nodeIdx,2)));
	hessAnalyticU(nodeIdx,2,1) = -2.0*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*(lambda2*sin(lambda2*nodePositions(nodeIdx,2))+lambda1*cos(lambda2*nodePositions(nodeIdx,2)));
	hessAnalyticU(nodeIdx,2,2) = -2.0*lambda2*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*(lambda2*cos(lambda2*nodePositions(nodeIdx,2))-lambda1*sin(lambda2*nodePositions(nodeIdx,2)))+2.0*lambda1*exp(nodePositions(nodeIdx,1))*exp(-lambda1*nodePositions(nodeIdx,2))*(lambda2*sin(lambda2*nodePositions(nodeIdx,2))+lambda1*cos(lambda2*nodePositions(nodeIdx,2)));
      endfor
    case STANDARD_LAPLACE_ANALYTIC_TYPE
      for nodeIdx = 1:numNodes
	analyticU(nodeIdx,1) = nodePositions(nodeIdx,1)*nodePositions(nodeIdx,1)+2*nodePositions(nodeIdx,1)*nodePositions(nodeIdx,2)-nodePositions(nodeIdx,2)*nodePositions(nodeIdx,2);
	gradAnalyticU(nodeIdx,1) = 2.0*nodePositions(nodeIdx,1)+2.0*nodePositions(nodeIdx,2);
	gradAnalyticU(nodeIdx,1) = 2.0*nodePositions(nodeIdx,1)-2.0*nodePositions(nodeIdx,2);
	hessAnalyticU(nodeIdx,1,1) = 0.0;
	hessAnalyticU(nodeIdx,1,2) = 0.0;
	hessAnalyticU(nodeIdx,2,1) = 0.0;
	hessAnalyticU(nodeIdx,2,2) = 0.0;
      endfor
    otherwise
      printf("ERROR: Invalid analytic type.\n");
      quit;
  endswitch
  
  if( diagnostics == DIAGNOSTICS_ON )
    analyticU
    #gradAnalyticU
    #hessAnalyticU
  endif
  
endfunction

function [ u, dofMap, freeDOFS, fixedDOFS ] = setBoundaryConditions( boundaryNodes, analyticU )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global ON_BOUNDARY_NODE;
  global DOF_FIXED;
  global DOF_FREE;
  global numNodes;

  dofIdx = 0;
  freeDOFIdx = 0;
  fixedDOFIdx = 0;
  for nodeIdx = 1:numNodes
    dofIdx = dofIdx + 1;
    if( boundaryNodes( nodeIdx ) == ON_BOUNDARY_NODE )
      fixedDOFIdx = fixedDOFIdx + 1;
      dofMap( dofIdx, 1 ) = DOF_FIXED;
      dofMap( dofIdx, 2 ) = fixedDOFIdx;
      fixedDOFS( fixedDOFIdx ) = dofIdx;
      u( nodeIdx, 1 ) = analyticU( nodeIdx, 1 );
    else
      freeDOFIdx = freeDOFIdx + 1;
      dofMap( dofIdx, 1 ) = DOF_FREE;
      dofMap( dofIdx, 2 ) = freeDOFIdx;
      freeDOFS( freeDOFIdx ) = dofIdx;
      u( nodeIdx, 1 ) = 0.0;
    endif
  endfor
  
  if( diagnostics == DIAGNOSTICS_ON )
    u
    dofMap
    freeDOFS
    fixedDOFS
  endif
  
endfunction

function [ A, b ] = reduceGlobalSystem( dofMap, freeDOFS, fixedDOFS, K, u, f )
  
  global DIAGNOSTICS_ON;
  global diagnostics;
  global DOF_FIXED;
  global DOF_FREE;
  global numDOFS;
  global numFreeDOFS;
  global numFixedDOFS;
  
  for rowDOFIdx = 1:numFreeDOFS
    globalRowNumber = freeDOFS( rowDOFIdx );
    # Form reduced A matrix
    for columnDOFIdx = 1:numFreeDOFS
      globalColumnNumber = freeDOFS( columnDOFIdx );
      A( rowDOFIdx, columnDOFIdx ) = K( globalRowNumber, globalColumnNumber );      
    endfor
    # Form RHS vector
    b( rowDOFIdx, 1 ) = f( globalRowNumber, 1 );
    for columnDOFIdx = 1:numFixedDOFS
      globalColumnNumber = fixedDOFS( columnDOFIdx );
      b( rowDOFIdx, 1 ) = b( rowDOFIdx, 1 ) - K( globalRowNumber, globalColumnNumber )*u( globalColumnNumber, 1 );
    endfor
  endfor

  if( diagnostics == DIAGNOSTICS_ON )
    A
    b
    u
    K*u
  endif
  
endfunction

function [ x ] = solveSystem( A, b );

  global DIAGNOSTICS_ON;
  global diagnostics;

  x = linsolve( A, b );

  if( diagnostics == DIAGNOSTICS_ON )
    x
  endif

endfunction

function [ u, f ] = updateDOFValues( dofMap, freeDOFS, fixedDOFS, x, K, u, f )
  
  global DIAGNOSTICS_ON;
  global diagnostics;
  global DOF_FIXED;
  global DOF_FREE;
  global numDOFS;
  global numFreeDOFS;
  global numFixedDOFS;
  #global u;
  #global f;
  global tempRHS=zeros(numDOFS,1)

  # Update u with solved values
  for dofIdx = 1:numFreeDOFS
    globalDOFIdx = freeDOFS( dofIdx );
    u( globalDOFIdx, 1 ) = x( dofIdx, 1 );
    f( globalDOFIdx, 1 ) = f( globalDOFIdx, 1 );
  endfor
  # Backsubstitute
  tempRHS = K*u;
  for dofIdx = 1:numFixedDOFS
    globalDOFIdx = fixedDOFS( dofIdx );
    u( globalDOFIdx, 1 ) = u( globalDOFIdx );
    f( globalDOFIdx, 1 ) = -1.0*tempRHS( globalDOFIdx, 1 );
  endfor
  
  if( diagnostics == DIAGNOSTICS_ON )
    u
    f
  endif
  
endfunction

function [ err, percentErr, rmsErr ] = computeErrors( u, analyticU )

  global DIAGNOSTICS_ON;
  global diagnostics;
  global numDOFS;
  
  sum = 0.0;
  printf("  row    value analytic      err    %%err \n");
  for dofIdx = 1:numDOFS
    err(dofIdx,1) = (u(dofIdx,1)-analyticU(dofIdx,1));
    if( abs(analyticU(dofIdx,1)) > 0.00000001 )
      percentErr(dofIdx,1) = 100.0*(u(dofIdx,1)-analyticU(dofIdx,1))/analyticU(dofIdx,1);
    else
      percentErr(dofIdx,1) = 0.0;
    endif
    sum = sum + err(dofIdx,1)*err(dofIdx,1);
    printf("%5d %8.5f %8.5f %8.5f %7.2f\n",dofIdx,u(dofIdx,1),analyticU(dofIdx,1),err(dofIdx,1),percentErr(dofIdx,1));
  endfor
  rmsErr = sqrt(sum/numDOFS)
  
  if( diagnostics == DIAGNOSTICS_ON )
  endif

endfunction

# Start of main code

global nodePositions = zeros(numNodes,2);
global boundaryNodes = zeros(numNodes,1);
global u = zeros(numNodes,1);
global analyticU = zeros(numNodes,1);
global gradAnalyticU = zeros(numNodes,2);
global hessAnalyticU = zeros(numNodes,2,2);
global K = zeros(numDOFS,numDOFS);
global f = zeros(numDOFS,1);
global dofMap = zeros(numDOFS,2);
global freeDOFS = zeros(numFreeDOFS,1);
global fixedDOFS = zeros(numFixedDOFS,1);
global A = zeros(numFreeDOFS,numFreeDOFS);
global x = zeros(numFreeDOFS,1);
global b = zeros(numFreeDOFS,1);
global uError = zeros(numNodes,1);
global uPerError = zeros(numNodes,1);
global uRMSError = 0.0;


[ nodePositions, boundaryNodes ] = computeGeometry( interpolationType );
[ analyticU, gradAnalyticU, hessAnalyticU ] = analytic( analyticType, nodePositions );
[ K, f ] = computeLinearMatrices( interpolationType );
[ u, dofMap, freeDOFS, fixedDOFS ] = setBoundaryConditions( boundaryNodes, analyticU );
[ A, b ] = reduceGlobalSystem( dofMap, freeDOFS, fixedDOFS, K, u, f );
[ x ] = solveSystem( A, b );
[ u, f ] = updateDOFValues( dofMap, freeDOFS, fixedDOFS, x, K, u, f );
[ err, percentErr, rmsErr ] = computeErrors( u, analyticU );
