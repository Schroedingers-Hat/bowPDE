  clear
  clf
  xPts = 101;                                                               % Number of x points. Odd
  L = 1.5;
  dx = L/xPts;
  dt   = 0.0000001;                                                              % Time step.
  m    = 1/100;                                                             % Mass density.
  mid = (xPts+1)/2;
  
  %coupled PDE vector
  cVec = [zeros(2*xPts,1)];

  % Set force function.
  q    = zeros(xPts, 1);                                                   % No force on most of it.
  q(5) = 1000;
  q(end-4) = 1000;

 q((xPts + 1)/2) = -2000;

  % Difference operator. 
  % Using second order coefficients for second derivative.

  d2Coeffs2 = [0 0 0 1 -2 1 0 0 0];  
  d1Coeffs1 = [0 0 0 -1 0 1 0 0 0];  
  d2Coeffs8 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560];
  dTwo = zeros(xPts,xPts);
  dOne = zeros(xPts,xPts);

  for count = -4:4
    % Add an offset diagonal matrix for each step to build banded matrix.
    dTwo = dTwo + ...
            d2Coeffs2(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
    dOne = dOne + ...
            d1Coeffs1(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
  end
  dTwo(1,1) = -1;
  dTwo(end,end) = -1;
  E = 1E9;
  I = [1:2/xPts:2];
  I = [I fliplr(I(1:end-1))];
  EI = diag(E.*I);
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 = -dt/m*EBOp - 10000 * dTwo;  %zeros(xPts,xPts);                                     % map v to v stepping backwards
  M2 = -EBOp/m;                                                            % map w to v
  M3 = eye(xPts);                                                          % map w to w
  M4 = zeros(xPts,xPts);                                                   % map v to w
  
  Tfwd = [M1 M2;...
          M3 M4];
  Tcrank = (eye(2*xPts) - Tfwd * dt)^-1;
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
  
cphi = zeros(xPts,1);
X = cphi;
Y = cphi;

  hold off;
  fin = 1000000;
  for count = 1:fin;
      
    dw = cVec(xPts+2:end) - cVec(xPts+1:end-1);
    dphi = atan(dw/dx);
    
    for count2 = 1:mid - 1
      cphi(mid - count2) = cphi(mid - count2 + 1) + dphi(mid - count2);
      cphi(mid + count2) = cphi(mid + count2 - 1) + dphi(mid + count2 -1);
      
      X(mid - count2) = X(mid - count2 + 1) + cos(cphi(mid - count2));
      X(mid + count2) = X(mid + count2 - 1) - cos(cphi(mid + count2));
      
      Y(mid - count2) = Y(mid - count2 + 1) - sin(cphi(mid - count2));
      Y(mid + count2) = Y(mid + count2 - 1) + sin(cphi(mid + count2));   
    end
    cphi = abs(cphi);
    
    F = [q/m*dt; zeros(xPts,1)].*cos([cphi; zeros(xPts,1)]);
    cVec = Tcrank * (cVec + F);
%     m1 = dt*Tfwd*cVec+F;
%     m2 = dt*Tfwd*cVec+m1/2+F;
%     m3 = dt*Tfwd*cVec+m2/2+F;
%     m4 = dt*Tfwd*cVec+m3+F;
%     
%     cVec = cVec + (1/6)*(m1+2*(m2+m3)+m4);
%     

    cVec = setBoundaries(cVec,xPts);
    
    if mod(count, fin/1000) == 0
        clc
        pc = count*100/fin
        hold off;
%         figure(1)
%         plot(cVec(xPts+1:2*xPts));
%         figure(2)
%         plot(F(1:xPts))
%         figure(3)
%         plot(cphi)
        figure(4)
        plot(X,Y)
        axis([-mid mid -15 50])
    end
  end
