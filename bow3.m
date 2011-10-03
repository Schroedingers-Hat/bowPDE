  clear
  clf
  xPts = 101;                                                               % Number of x points. Odd
  dt   = 0.001;                                                              % Time step.
  m    = 1/10;                                                             % Mass density.
  
  %coupled PDE vector
  cVec = [zeros(2*xPts,1)];

  % Set force function.
  q    = zeros(xPts, 1);                                                   % No force on most of it.
  q(5) = 0.001;
  q(end-4) = 0.001;
%   q((xPts + 1)/2) = -0.002;

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
  I = [1:2/xPts:2];
  I = [I fliplr(I(1:end-1))];
  EI = diag(E.*I);
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 = -dt/m*EBOp;  %zeros(xPts,xPts);                                                      % map v to v stepping backwards
  M2 = -EBOp/m;                                                               % map w to v
  M3 = eye(xPts);                                                           % map w to w
  M4 = zeros(xPts,xPts);                                                    % map v to w
  
  Tfwd = [M1 M2;...
          M3 M4];
 
%   Boundary conditions.
%   Tfwd((xPts + 1) / 2,:) = 0;
%   Tfwd((xPts + 1) / 2,(xPts + 1) / 2) = 1;      %Avoid singular matrix
%   Tfwd((3*xPts + 1) / 2,:) = 0;
%   Tfwd((3*xPts + 1) / 2,(3*xPts + 1) / 2) = 1;  %Avoid singular matrix
%   Tfwd(1:4,:) = 0;
%   Tfwd(1:4,5) = 1;                             
%   Tfwd(xPts-3:xPts,:) = 0;
%   Tfwd(xPts-3:xPts,xPts-4) = 1;
%   
%   bw = 4;
%   for k = 1:bw
%     Tfwd(xPts + bw + 1 - k,:) = 0;
%     Tfwd(xPts + bw + 1 - k,xPts + bw + 2) = 1-k;
%     Tfwd(xPts + bw + 1 - k,xPts + bw + 1) = k;
%     Tfwd(end - bw + k,:) = 0;
%     Tfwd(end - bw + k,end - bw - 1) = 1-k;
%     Tfwd(end - bw + k,end - bw) = k;
%   end
  
  % Do similar things for the straight ends.
  % Could use lower order finite difference eqs so only the last few pts


  % Won't bother combining things into a single matrix yet.
  hold off;
  fin = 100000;
  for count = 1:fin;
      
    %cVec = Tfwd * cVec;
    F = [q/m*dt; zeros(xPts,1)];
    
    m1 = dt*Tfwd*cVec+F;
    m2 = dt*Tfwd*cVec+m1/2+F;
    m3 = dt*Tfwd*cVec+m2/2+F;
    m4 = dt*Tfwd*cVec+m3+F;
    
    cVec = cVec + (1/6)*(m1+2*(m2+m3)+m4);
    
    cVec = setBoundaries(cVec,xPts);
    
    if mod(count, fin/100) == 0
        clc
        pc = count*100/fin
        hold on;
        plot(cVec(xPts+1:2*xPts));
    end
  end
