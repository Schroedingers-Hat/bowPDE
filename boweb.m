  clear
  clf
  xPts = 305;                                                               % Number of x points. Odd
  L = 2;
  dx = L/xPts;
  dt   = 0.001;                                                              % Time step.
  rho    = 670;                                                             % Mass density.
  mid = (xPts+1)/2;
  k = 5E6;                                                                  %springyness of string
  Fx = [];
  Fy = [];
  m_arrow = 0.3;
  
  E = 1.5E10;
  a = 1E-2*[2*ones(1,xPts/5) 2.5*ones(1,xPts/5) 3.2*ones(1,xPts/5) 2.5*ones(1,xPts/5) 2*ones(1,xPts/5)];                                                   %width
  b = 1E-2*[2*ones(1,xPts/5) 3.4*ones(1,xPts/5) 4*ones(1,xPts/5) 3.4*ones(1,xPts/5) 2*ones(1,xPts/5)];                                                   %width;                                                  %thickness

  md = diag(a.*b.*rho.*dx);
  lmd = [md^-1 zeros(xPts);zeros(xPts) zeros(xPts)];
  
  I = (1/3)*a.*b.^3;                                                       %moment of area
  EI = diag(E.*I);
  
  arr = [-1:2/(xPts-1):1];
  deflection = zeros(xPts,1);
   deflection = 0.*sinc(arr*3)';         %dw/dx

  %initial conditions
  cVec = zeros(2*xPts,1);
  
  % Set force function.
  q    = zeros(xPts,1); % No force on most of it.
  fk = 1;
  q(5) = 1*fk;%*dx^4;
  q(end-4) = 1*fk;%*dx^4;
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
  dTwo(1,1) = -1;
  dTwo(end,:) = 0;
  dTwo(end,end) = -1;
  
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 = -1*dt*md^(-1)*EBOp;                                                % map v to v stepping backwards
  M2 = -md^(-1)*EBOp;                                                          % map w to v
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

  while Y(5) <= 0.15
     
coord_transform
    
    if fap == 1
          woot = 0.001;
          sl = X(end - 4)*0.9;
          sy0 = Y(5);
          sy = sy0;
    end
    
    
    sdx = X(end-4);
    sdy = sy - Y(5);
    sdl = sqrt(sdx^2+sdy^2)-sl;
    sa = atan((sy-Y(5))/sdx);

    q = q+0.1*q0;
    F = [q;zeros(xPts,1)];
   
    cVec = Tcrank * (cVec + lmd*F);                                             %CRANK


        if mod(fap, 10) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        force_on_arm = F(5)
        applied_force = -q(mid)*dt
        hold off;
        figure(1)
        plot(Y,X)
        hold on
        axis equal
        axis([-0.5 1 -1.2 1.2])
        end
    fap = fap+1;
  end
  cVec(1:xPts) = 0;
        figure(2)
        plot(Y,X)
        hold on
        axis equal
        axis([-0.2 0.5 -1.2 1.2])
%CALCULATION LOOP----------------------------------------------------------
count = 1;
hold off;
Fdx = [];
Fdy = [];
woot = 0.001;
sl = X(end - 4)*0.95;
sy00 = Y(5);
sy0 = Y(5);
sy = sy0;
go = 0;
    Fapplied = 0;
    sv = [0] 
    SS = [sv sy];
  while (go <=12/dt)
    
coord_transform

    
    if (sy >= Y(5)*1.01 && go <= 1/dt)
    sy = sy-1*dt*sy0;
    sy0 = sy;
    elseif (go > 1000 & go < 10/dt & sy < 0.7)
    sy = sy+0.5*dt*sy00;
    Fx = [Fx sy-sy0];
    Fy = [Fy -Fapplied];
    elseif (go > 10/dt & sy>sy0)
        go
        sy = sy + sv(end)*dt;
        sv = [sv sv(end)+Fapplied/m_arrow*dt];
        Fdx = [Fdx sy-sy0];
        Fdy = [Fdy -Fapplied];
    end
    go = go+1;
    sdx = X(end-4);
    sdy = sy - Y(5);
    sdl = max([sqrt(sdx^2+sdy^2)-sl 0]);
    sa = atan((sy-Y(5))/sdx);
    
    F = [max([k*sdl*sin(sa-cphi(5))*dt; woot]).*q0 ;zeros(xPts,1)];
%    F(mid) = -2*F(5);
    F(mid) = -2*k*sdl*max([sin(sa-cphi(5)) 0])*dt;
    Fapplied = -2*k*sdl*max([sin(sa) 0])*dt;
    cVec = Tcrank * (cVec + lmd*F);                                             %CRANK

    
    if mod(count, fin/100) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        force_on_arm = F(5)
        applied_force = -Fapplied
        stretch = sdl/sl
        sl
        hold off;

        figure(1)
        plot([Y(5) sy Y(end-4)],[X(5) 0 X(end-4)])
        hold on
        plot(Y(1:end/5),X(1:end/5),'k')
        plot(Y(end/5:2*end/5),X(1*end/5:2*end/5),'r')
        plot(Y(2*end/5:3*end/5),X(2*end/5:3*end/5),'k')
        plot(Y(3*end/5:4*end/5),X(3*end/5:4*end/5),'r')
        plot(Y(4*end/5:end),X(4*end/5:end),'k')
        axis equal
        axis([-0.2 1 -1.2 1.2])
    end
    count = count+1;
  end

  figure(3)
  plot(Fx, Fy)
hold on
plot([Fx(1),Fx(end)], [Fy(1),Fy(end)],'r')
  figure(69)
  plot(Fdx, Fdy)