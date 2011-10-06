  clear
  clf
  xPts = 405;                                                               % Number of x points. Odd
  L = 2;
  dx = L/xPts;
  dt   = 0.001;                                                              % Time step.
  m    = 1/xPts;                                                             % Mass density.
  mid = (xPts+1)/2;
  k = 1E6;                                                                  %springyness of string
  Fx = [];
  Fy = [];
  
  E = 1E9;
  a = 1E-2*[0.5*ones(1,xPts/5) 1*ones(1,xPts/5) 2*ones(1,xPts/5) 1*ones(1,xPts/5) 0.5*ones(1,xPts/5)];                                                   %width
  b = 1E-2*[0.5*ones(1,xPts/5) 0.5*ones(1,xPts/5) 1.5*ones(1,xPts/5) 0.5*ones(1,xPts/5) 0.5*ones(1,xPts/5)];                                                   %width;                                                  %thickness
  I = (1/12)*a.*b.^3;                                                       %moment of area
  EI = diag(E.*I);
  
  arr = [-1:2/(xPts-1):1];
  deflection = zeros(xPts,1);
  deflection(mid-floor(xPts/5):mid+floor(xPts/5)) = -0.2*arr(mid-floor(xPts/5):mid+floor(xPts/5)).^2';         %dw/dx
  deflection(1:mid-floor(xPts/5)-1) = deflection(mid-floor(xPts/5)) - ...
      diff(deflection(mid-floor(xPts/5)-1:mid-floor(xPts/5)))*arr(1:mid-floor(xPts/5)-1);         %dw/dx
  deflection(mid+floor(xPts/5)+1:end) = deflection(mid+floor(xPts/5)) - ...
      diff(deflection(mid+floor(xPts/5):mid+floor(xPts/5)+1))*arr(mid+floor(xPts/5)+1:end);         %dw/dx
  %initial conditions
  cVec = zeros(2*xPts,1);
  
  % Set force function.
  q    = zeros(xPts,1); % No force on most of it.
  fk = 1;
  q(3) = 1*fk;%*dx^4;
  q(end-3) = 1*fk;%*dx^4;
  q(mid) = -2*fk;%*dx^4;
  q0 = q;

  % Difference operator. 
  % Using second order coefficients for second derivative.

  d2Coeffs2 = [0 0 0 1 -2 1 0 0 0];  
  d1Coeffs = [0 0 0 -1/2 0 1/2 0 0 0];  
  dTwo = zeros(xPts,xPts);
  dOne = zeros(xPts,xPts);
  for count = -4:4
    % Add an offset diagonal matrix for each step to build banded matrix.
    dTwo = dTwo + ...
            (1/dx^2)*d2Coeffs2(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
    dOne = dOne + ...
            (1/dx)*d1Coeffs(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
  end
  dOne(1,1:3) = [-1 1 0]/dx;
  dOne(end,end-2:end) = [0 -1 1]/dx;
  dTwo(1,:) = 0;
  dTwo(1,1:3) = [1 -2 1];
  dTwo(end,:) = 0;
  dTwo(end,end-2:end) = [1 -2 1];
  
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 = -1000*dt/m*EBOp;                                                % map v to v stepping backwards
  M2 = -EBOp/(m);                                                          % map w to v
  M3 = eye(xPts);                                                          % map w to w
  M4 = zeros(xPts,xPts);                                                   % map v to w
  
  
  Tfwd = [M1 M2;...
          M3 M4];
  Tcrank = (eye(2*xPts) - Tfwd * dt)^-1;
  
cphi = zeros(xPts,1);
X = cphi;
Y = cphi;

%STRINGING LOOP------------------------------------------------------------
  fin = 10000;
  fap = 1;

  while Y(5) <= 0.2
     
coord_transform
    
    if fap == 1
          woot = 0.001;
          sl = X(end - 4)*0.9;
          sy0 = Y(5);
          sy = sy0;
    end
    
    if sy <= 1.4
    sy = sy+0.002*sy0;   
    end
    
    sdx = X(end-4);
    sdy = sy - Y(5);
    sdl = sqrt(sdx^2+sdy^2)-sl;
    sa = atan((sy-Y(5))/sdx);

    q = q+0.1*q0;
    F = [q;zeros(xPts,1)];
   
    cVec = Tcrank * (cVec + F);                                             %CRANK

    
        if mod(fap, 10) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        force_on_arm = F(5)
        applied_force = -q(mid)/m*dt
        sdx = sdx
        sa = sa
        hold off;
        figure(1)
        plot([Y(5)+sdy; Y(5); Y; Y(end-4); Y(5)+sdy],[0; X(5); X; X(end-4); 0]')
        hold on
        axis equal
        axis([-0.5 1 -1 1])
        end
    fap = fap+1;
  end
  cVec(1:xPts) = 0;
        figure(376)
        plot([Y(5); Y(5); Y; Y(end-4); Y(5)],[0; X(5); X; X(end-4); 0]')
        hold on
        axis equal
        axis([-0.2 2 -1 1])
%CALCULATION LOOP----------------------------------------------------------
  hold off;
  for count = 1:fin;
     
coord_transform
    
    if count == 1
          woot = 0.001;
          sl = X(end - 4);
          sy0 = Y(5);
          sy = sy0;
    end
    
    if sy <= 1.4 && count >= fin/10
    sy = sy+dt*sy0;   
    end
    
    sdx = X(end-4);
    sdy = sy - Y(5);
    sdl = sqrt(sdx^2+sdy^2)-sl;
    sa = atan((sy-Y(5))/sdx);
    
    F = [max([k*sdl*sin(sa-cphi(5))*dt; 0]).*q ;zeros(xPts,1)];
    F(mid) = -2*k*sdl*sin(sa)*dt;
    
    cVec = Tcrank * (cVec + F);                                             %CRANK

    
    Fx(count) = sy;
    Fy(count) = -F(mid);
    
    if mod(count, fin/100) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        force_on_arm = F(5)
        applied_force = -F(mid)
        sdl
        hold off;

        figure(1)
        plot([sy; Y(5); Y; Y(end-4); sy],[0; X(5); X; X(end-4); 0]')
        axis equal
        axis([-0.5 2 -1 1])
    end
  end

  figure
  plot(Fx, Fy)