  clear
  clf
  xPts = 105;                                                               % Number of x points. Odd
  L = 2;
  dx = L/xPts;
  dt   = 0.0001;                                                              % Time step.
  m    = 1/xPts;                                                             % Mass density.
  mid = (xPts+1)/2;
  k = 1E6;                                                                  %springyness of string
  
  E = 5E9*[1*ones(1,xPts/5) 1*ones(1,xPts/5) 1*ones(1,xPts/5) 1*ones(1,xPts/5) 1*ones(1,xPts/5)];
  a = ones(1,xPts)*3E-2;                                                    %width
  b = ones(1,xPts)*0.5E-2;                                                  %thickness
  I = (1/12)*a.*b.^3;                                                       %moment of area
  EI = diag(E.*I);
  deflection = -000*0.000154311829612253*ones(xPts,1);         %dw/dx

  %deflection(mid-100:mid+100) = -00*[0.000154311829612253;0.000151350390876264;0.000148415486115593;0.000145507139697444;0.000142625289098539;0.000139770264119298;0.000136941944480411;0.000134140397811805;0.000131365744732586;0.000128618032825272;0.000125897432456692;0.000123203842704522;0.000120537440671360;0.000117898225286510;0.000115286275826947;0.000112701740391245;0.000110144507447008;0.000107614808576741;0.000105112658598623;0.000102638130509692;0.000100191302091868;9.77722191128426e-05;9.53808964385729e-05;9.30174481231315e-05;9.06820615646065e-05;8.83745854724635e-05;8.60952700526999e-05;8.38439591858228e-05;8.16208983437442e-05;7.94261800740458e-05;7.72597612310136e-05;7.51217923998213e-05;7.30123210107600e-05;7.09312539870368e-05;6.88788986897588e-05;6.68553100002336e-05;6.48602981508133e-05;6.28941512460368e-05;6.09569212821371e-05;5.90486080268472e-05;5.71692961350332e-05;5.53190477279478e-05;5.34979234469818e-05;5.17060653055418e-05;4.99433974416159e-05;4.82100559035391e-05;4.65062162866911e-05;4.48317775242905e-05;4.31869549805985e-05;4.15717092840138e-05;3.99860560365259e-05;3.84302126274932e-05;3.69041214208592e-05;3.54078744393225e-05;3.39416839558766e-05;3.25054099401998e-05;3.10991421562604e-05;2.97230949939253e-05;2.83773561876881e-05;2.70616927387921e-05;2.57764142670203e-05;2.45215363404740e-05;2.32971930083153e-05;2.21033988239082e-05;2.09401088140099e-05;1.98075270254398e-05;1.87056122252137e-05;1.76345830532813e-05;1.65943450197397e-05;1.55850855539228e-05;1.46067974903564e-05;1.36595807570760e-05;1.27435388067614e-05;1.18585682572438e-05;1.10048553748246e-05;1.01826818814754e-05;9.39179139157509e-06;8.63228830347989e-06;7.90438135485305e-06;7.20805751124317e-06;6.54334080573339e-06;5.91038839766147e-06;5.30918218732123e-06;4.73994044808096e-06;4.20255119766201e-06;3.69713930686081e-06;3.22390899571897e-06;2.78258469668642e-06;2.37363011529491e-06;1.99679845764130e-06;1.65235813272959e-06;1.34030412524977e-06;1.06069419323439e-06;8.13587596908364e-07;5.99021699200739e-07;4.17184642597914e-07;2.68132705995560e-07;1.51742072120153e-07;6.82428450962485e-08;1.76800765473787e-08;0;1.54269251260599e-08;6.39692816081828e-08;1.45311175039842e-07;2.59507938745511e-07;4.06433679463086e-07;5.86219309902757e-07;7.98455009395647e-07;1.04343556783390e-06;1.32088845552159e-06;1.63079539093618e-06;1.97315367215801e-06;2.34774166484022e-06;2.75453854106646e-06;3.19364000013485e-06;3.66483417291850e-06;4.16807474236621e-06;4.70329046616366e-06;5.27050624554910e-06;5.86937793586805e-06;6.50025643975950e-06;7.16272725150731e-06;7.85694208933386e-06;8.58265385862717e-06;9.34004508366230e-06;1.01287810431404e-05;1.09489536347158e-05;1.18003853339062e-05;1.26831164520234e-05;1.35970692838977e-05;1.45421466249100e-05;1.55182442094690e-05;1.65253874893998e-05;1.75634678519198e-05;1.86323038711294e-05;1.97320783370865e-05;2.08625266154302e-05;2.20236141828255e-05;2.32152731037684e-05;2.44375245703211e-05;2.56902201977524e-05;2.69733385008382e-05;2.82867779666510e-05;2.96304420471589e-05;3.10043372564382e-05;3.24083420975188e-05;3.38425133444917e-05;3.53065457890846e-05;3.68005715186669e-05;3.83245389358488e-05;3.98782653882510e-05;4.14616952164644e-05;4.30748188780442e-05;4.47175650478249e-05;4.63898467993801e-05;4.80915319034599e-05;4.98226971677455e-05;5.15831512232700e-05;5.33728759952781e-05;5.51918321522111e-05;5.70398932535181e-05;5.89170846351397e-05;6.08231974387505e-05;6.27583343299031e-05;6.47222792528885e-05;6.67151027174107e-05;6.87366342110330e-05;7.07868142156144e-05;7.28656406545752e-05;7.49730151333333e-05;7.71088514536503e-05;7.92731023028393e-05;8.14657486020363e-05;8.36865503218298e-05;8.59356086301121e-05;8.82128718954112e-05;9.05181140965237e-05;9.28514356830927e-05;9.52127195198000e-05;9.76018480893283e-05;0.000100018833389875;0.000102463466950951;0.000104935888582002;0.000107435853506894;0.000109963427809404;0.000112518416931519;0.000115100878417533;0.000117710619245428;0.000120347656556815;0.000123011958360559;0.000125703392813205;0.000128421904837227;0.000131167359042717;0.000133939843777353;0.000136739170060611;0.000139565385326211;0.000142418283302892;0.000145297878316896;0.000148204149921329;0.000151136894680200;0.000154096147220371;];
 
  
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
  dTwo = zeros(xPts,xPts);

  for count = -4:4
    % Add an offset diagonal matrix for each step to build banded matrix.
    dTwo = dTwo + ...
            (1/dx^2)*d2Coeffs2(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
  end
  dTwo(1,:) = 0;
  dTwo(1,1) = 0;
  dTwo(end,:) = 0;
  dTwo(end,end) = 0;
  
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 = -10*dt/m*EBOp;                                                % map v to v stepping backwards
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
  fin = 1000;
  fap = 1;
  while Y(5) <= 0.2;
     
    dw1 = diff(cVec(xPts+1:end) + deflection);
    dw2 = diff(flipud(cVec(xPts+1:end) + deflection));
    phi1 = (atan(dw1/dx));
    phi2 = (atan(dw2/dx));
    phi = (phi1+phi2)/2;
    dphi = [0;diff(phi);0];
    
    for count2 = 1:mid - 1
      cphi(mid - count2) = cphi(mid - count2 + 1) + dphi(mid - count2);
      cphi(mid + count2) = cphi(mid + count2 - 1) + dphi(mid + count2 -1);
      
      X(mid - count2) = X(mid - count2 + 1) - dx*cos(cphi(mid - count2));
      X(mid + count2) = X(mid + count2 - 1) + dx*cos(cphi(mid + count2));
      
      Y(mid - count2) = Y(mid - count2 + 1) + dx*sin(cphi(mid - count2));
      Y(mid + count2) = Y(mid + count2 - 1) + dx*sin(cphi(mid + count2));   
    end
    
    if fap == 1
          woot = 0.01;
          sl = X(end - 4);
          sy0 = Y(5);
    end
    
    sdx = X(end - 4);        %string x component
    sa = acos(sdx/sl);              %angle of string
    sdy = sl*sin(sa);             %string y component

    q = q+0.1*q0;
   F = abs([q/m*dt; zeros(xPts,1)]*sin((-cphi(5)+sa)/(2*max([(sin(sa)) woot]))));
   F(mid) = q(mid)/m*dt;
   
    cVec = Tcrank * (cVec + F);                                             %CRANK
%     cVec(xPts+1:2*xPts) = cVec(xPts+1:2*xPts) - cVec(mid+xPts);             %move back
%     cVec(1:xPts) = cVec(1:xPts) - cVec(mid); 
    
        if mod(count, fin/100) == 0
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
    
    fap = fap + 0.000001
  end
  cVec(1:xPts) = 0;
%         figure(376)
%         plot([Y(5); Y(5); Y; Y(end-4); Y(5)],[0; X(5); X; X(end-4); 0]')
%         hold on
%         axis equal
%         axis([-0.2 sdy+0.2 -1 1])
%CALCULATION LOOP----------------------------------------------------------
  hold off;
  for count = 1:fin;
     
    dw1 = diff(cVec(xPts+1:end) + deflection);
    dw2 = diff(flipud(cVec(xPts+1:end) + deflection));
    phi1 = (atan(dw1/dx));
    phi2 = (atan(dw2/dx));
    phi = (phi1+phi2)/2;
    dphi = [0;diff(phi);0];
    
    for count2 = 1:mid - 1
      cphi(mid - count2) = cphi(mid - count2 + 1) + dphi(mid - count2);
      cphi(mid + count2) = cphi(mid + count2 - 1) + dphi(mid + count2);
      
      X(mid - count2) = X(mid - count2 + 1) - dx*cos(cphi(mid - count2));
      X(mid + count2) = X(mid + count2 - 1) + dx*cos(cphi(mid + count2));
      
      Y(mid - count2) = Y(mid - count2 + 1) + dx*sin(cphi(mid - count2));
      Y(mid + count2) = Y(mid + count2 - 1) + dx*sin(cphi(mid + count2));   
    end
    
    if count == 1
          woot = 0.001;
          sl = X(end - 4);
          sy0 = Y(5);
    end

    sy = count/fin+sy0;                                        %draw distance
    stretch = sqrt((X(mid)-X(5)).^2+(sy-Y(5)).^2)-sl
    
    sdx = X(end - 4);        %string x component
    sa = atan((sy-sy0)/sdx);              %angle of string
    sdy = sl*sin(sa);             %string y component
    
    F = ([k*stretch*max(0,sin(sa-cphi(5)))/2*q0;zeros(xPts,1)]);
    F(mid) = -(k*stretch*max(0,sin(cphi(5))));
   
    cVec = Tcrank * (cVec + F);                                             %CRANK
    cVec(xPts+1:2*xPts) = cVec(xPts+1:2*xPts) - cVec(mid+xPts);             %move back
    cVec(1:xPts) = cVec(1:xPts) - cVec(mid);                                %decelerate
    
    %cVec = setBoundaries(cVec,xPts);
    
    if mod(count, fin/100000) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        force_on_arm = F(5)
        applied_force = -F(mid)
        sdx = sdx
        sa = sa
        hold off;
%        figure(2)
%         plot(cVec(xPts+1:2*xPts));
%        figure(72)
%         plot(dw);
%        axis([0 xPts -1E-3 1E-2])
%         figure(2)
%         plot(F(1:xPts))
%         figure(3)
%         plot(dphi)
%         axis([0 xPts -0.1 0.1])
         figure(4)
        plot(sa)
        axis([0 xPts -1 10])
        figure(1)
        plot([sy; Y(5); Y; Y(end-4); sy],[0; X(5); X; X(end-4); 0]')
        axis equal
        axis([-0.5 1 -1 1])
    end
  end
