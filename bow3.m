  clear
  clf
  xPts = 101;                                                               % Number of x points. Odd
  dt   = 0.01;                                                              % Time step.
  m    = 1/10;                                                             % Mass density.
  
  %coupled PDE vector
  cVec = [zeros(2*xPts,1); 1];

  % Set force function.
  q    = zeros(xPts, 1);                                                   % No force on most of it.
  q(5) = 0.001;
  q(end-4) = 0.001;
%  q((xPts + 1)/2) = 0.01;

  % Difference operator. 
  % Using second order coefficients for second derivative.

  d2Coeffs2 = [0 0 0 1 -2 1 0 0 0];  
  d2Coeffs8 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560];
  dTwo = zeros(xPts,xPts);

  for count = -4:4
    % Add an offset diagonal matrix for each step to build banded matrix.
    dTwo = dTwo + ...
            d2Coeffs2(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
  end
  dTwo(1,1) = -1;
  dTwo(end,end) = -1;

  E = 1;
  I = [1:2/xPts:2 2-2/xPts:-2/xPts:1].^4;
  EI = diag(E.*I);
  
  %EBOp = dFour;
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tbck U(x,t+1) = U(x,t)
%  M1 = -dt^2/m*EBOp+eye(xPts);                                            % map v to v
  
  M1 = eye(xPts);                                                          % map v to v
  M2 = EBOp/m*dt;                                                          % map w to v
  M3 = -(q/m*dt);                                                          % map 1 to v
  M4 = -eye(xPts)*dt;                                                      % map w to w
  M5 = eye(xPts);                                                          % map v to w
  M6 = zeros(xPts,1);                                                      % map 1 to w
  M7 = [zeros(1,2*xPts), 1];                                               % map 1 to 1
  
  Tbck = [M1 M2 M3;...
          M4 M5 M6;...
             M7];
  % Boundary conditions.
  Tbck((xPts + 1) / 2,:) = 0;
  Tbck((xPts + 1) / 2,(xPts + 1) / 2) = 1;
  Tbck((3*xPts + 1) / 2,:) = 0;
  Tbck((3*xPts + 1) / 2,(3*xPts + 1) / 2) = 1;
  % Do similar things for the straight ends.
  % Could use lower order finite difference eqs so only the last few pts

%[v,w](x,t+1) = (I - OPERATOR*dt)^(-1)[v,w](x,t)

  % Invert TOp to get forward time operator
  Tfwd = ((((Tbck^(-1)))^10)^200);
  % Won't bother combining things into a single matrix yet.
  hold off;
  for count = 1:10;
    cVec = Tfwd * cVec;
    cVec = setBoundaries(cVec,xPts);
    hold on;
    plot(cVec(xPts+1:2*xPts));
  end
