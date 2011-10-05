  clear
  clf
  xPts = 101;                                                               % Number of x points. Odd
  dt   = 0.0001;                                                              % Time step.
  m    = 1/10;                                                             % Mass density.
  
  %coupled PDE vector
  cVec = [zeros(4*xPts,1)];

  % Set force function.
  q    = zeros(xPts, 1);                                                   % No force on most of it.
%   q(5) = 0.001;
%   q(end-4) = 0.001;
  q((end+1)/2) = 0.001;
  
  % Difference operators. 

  d2Coeffs = [0 1/90 -3/20 3/2 -49/18 3/2 -3/20 1/90 0];
  d1Coeffs = [0 -1/60 3/20 -3/4 0 3/4 -3/20 1/60 0];
  dTwo = zeros(xPts,xPts);
  dOne = dTwo;
 
  
  for count = -4:4
    % Add an offset diagonal matrix for each step to build banded matrix.
    dTwo = dTwo + ...
            d2Coeffs(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
    dOne = dOne + ...
            d1Coeffs(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
  end
  dTwo(1,1) = -1;
  dTwo(end,end) = -1;
  
  % Set Constants
  a = [1:2/xPts:2];
  a = [a fliplr(a(1:end-1))];
  b = [1:2/xPts:2];
  b = [b fliplr(b(1:end-1))];
  E = 1E9;
  I = diag(a.*b.^3);
  p = 5;
  k = 1;
  G = 1;
  A = diag(a.*b);
  kAG = k*G*A;
  EI = E*I;
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1) [v w]
%   M1 = -dt/m*EBOp; %zeros(xPts,xPts);                                      % map v to v stepping backwards
%   M2 = -EBOp/m;                                                            % map w to v
%   M3 = eye(xPts);                                                          % map v to w
%   M4 = zeros(xPts,xPts);                                                   % map w to w
%   
%   Tfwd = [M1 M2;...
%           M3 M4];
 
  %Make HUGE MATRIX for timoshenko equation [v w dphi phi]
  %{
  v = dw/dt
  d/dt[v w dphi phi] = [(1/pA)d/dx(kAG(dw/dx-phi) v]
  %}
  M1 = zeros(xPts,xPts);                                    
  M2 = (1./p)*A^(-1)*dOne*(kAG*(dOne));
  M3 = zeros(xPts,xPts);
  M4 = (1./p)*A^(-1)*dOne*(-kAG);
  
  M5 = eye(xPts);
  M6 = zeros(xPts,xPts);
  M7 = zeros(xPts,xPts);
  M8 = zeros(xPts,xPts);
  
  M9 = zeros(xPts,xPts);
  M10 = 1/p*I^(-1)*kAG*(dOne);
  M11 = zeros(xPts,xPts);
  M12 = 1/p*I^(-1)*(dOne*EI*dOne-kAG);
  
  M13 = zeros(xPts,xPts);
  M14 = zeros(xPts,xPts);
  M15 = eye(xPts);
  M16 = zeros(xPts,xPts);
  
  Tfwd = [M1  M2  M3  M4;...
          M5  M6  M7  M8;...
          M9  M10 M11 M12;...
          M13 M14 M15 M16];

  hold off;
  fin = 10000;
  %F is force for [v w dphi phi]
  for count = 1:fin;
    F = [q./(p.*a'.*b')*dt; zeros(3*xPts,1)];
    
    m1 = dt*Tfwd*cVec+F;
    m2 = dt*Tfwd*cVec+m1/2+F;
    m3 = dt*Tfwd*cVec+m2/2+F;
    m4 = dt*Tfwd*cVec+m3+F;
    
    cVec = cVec + (1/6)*(m1+2*(m2+m3)+m4);
    
    %cVec = setBoundaries(cVec,xPts);
    
    if mod(count, fin/100) == 0
        clc
        pc = count*100/fin
        hold on;
        plot(cVec(xPts+1:2*xPts));
    end
  end
