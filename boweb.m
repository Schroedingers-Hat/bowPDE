  clear
  clf
  xPts = 1001;                                                               % Number of x points. Odd
  L = 1.5;
  dx = L/xPts;
  dt   = 0.0001;                                                              % Time step.
  m    = 1/1000;                                                             % Mass density.
  mid = (xPts+1)/2;
  
  E = 1E9*[ones(1,xPts/5) 1*ones(1,xPts/5) ones(1,xPts/5+1) 1*ones(1,xPts/5) ones(1,xPts/5)];
  a = ones(1,xPts);
  b = ones(1,xPts);
  I = (1/12)*a.*b.^3;
  EI = diag(E.*I);
  deflection = 5*0.000154311829612253*ones(xPts,1);         %dw/dx

  deflection(mid-100:mid+100) = 5*[0.000154311829612253;0.000151350390876264;0.000148415486115593;0.000145507139697444;0.000142625289098539;0.000139770264119298;0.000136941944480411;0.000134140397811805;0.000131365744732586;0.000128618032825272;0.000125897432456692;0.000123203842704522;0.000120537440671360;0.000117898225286510;0.000115286275826947;0.000112701740391245;0.000110144507447008;0.000107614808576741;0.000105112658598623;0.000102638130509692;0.000100191302091868;9.77722191128426e-05;9.53808964385729e-05;9.30174481231315e-05;9.06820615646065e-05;8.83745854724635e-05;8.60952700526999e-05;8.38439591858228e-05;8.16208983437442e-05;7.94261800740458e-05;7.72597612310136e-05;7.51217923998213e-05;7.30123210107600e-05;7.09312539870368e-05;6.88788986897588e-05;6.68553100002336e-05;6.48602981508133e-05;6.28941512460368e-05;6.09569212821371e-05;5.90486080268472e-05;5.71692961350332e-05;5.53190477279478e-05;5.34979234469818e-05;5.17060653055418e-05;4.99433974416159e-05;4.82100559035391e-05;4.65062162866911e-05;4.48317775242905e-05;4.31869549805985e-05;4.15717092840138e-05;3.99860560365259e-05;3.84302126274932e-05;3.69041214208592e-05;3.54078744393225e-05;3.39416839558766e-05;3.25054099401998e-05;3.10991421562604e-05;2.97230949939253e-05;2.83773561876881e-05;2.70616927387921e-05;2.57764142670203e-05;2.45215363404740e-05;2.32971930083153e-05;2.21033988239082e-05;2.09401088140099e-05;1.98075270254398e-05;1.87056122252137e-05;1.76345830532813e-05;1.65943450197397e-05;1.55850855539228e-05;1.46067974903564e-05;1.36595807570760e-05;1.27435388067614e-05;1.18585682572438e-05;1.10048553748246e-05;1.01826818814754e-05;9.39179139157509e-06;8.63228830347989e-06;7.90438135485305e-06;7.20805751124317e-06;6.54334080573339e-06;5.91038839766147e-06;5.30918218732123e-06;4.73994044808096e-06;4.20255119766201e-06;3.69713930686081e-06;3.22390899571897e-06;2.78258469668642e-06;2.37363011529491e-06;1.99679845764130e-06;1.65235813272959e-06;1.34030412524977e-06;1.06069419323439e-06;8.13587596908364e-07;5.99021699200739e-07;4.17184642597914e-07;2.68132705995560e-07;1.51742072120153e-07;6.82428450962485e-08;1.76800765473787e-08;0;1.54269251260599e-08;6.39692816081828e-08;1.45311175039842e-07;2.59507938745511e-07;4.06433679463086e-07;5.86219309902757e-07;7.98455009395647e-07;1.04343556783390e-06;1.32088845552159e-06;1.63079539093618e-06;1.97315367215801e-06;2.34774166484022e-06;2.75453854106646e-06;3.19364000013485e-06;3.66483417291850e-06;4.16807474236621e-06;4.70329046616366e-06;5.27050624554910e-06;5.86937793586805e-06;6.50025643975950e-06;7.16272725150731e-06;7.85694208933386e-06;8.58265385862717e-06;9.34004508366230e-06;1.01287810431404e-05;1.09489536347158e-05;1.18003853339062e-05;1.26831164520234e-05;1.35970692838977e-05;1.45421466249100e-05;1.55182442094690e-05;1.65253874893998e-05;1.75634678519198e-05;1.86323038711294e-05;1.97320783370865e-05;2.08625266154302e-05;2.20236141828255e-05;2.32152731037684e-05;2.44375245703211e-05;2.56902201977524e-05;2.69733385008382e-05;2.82867779666510e-05;2.96304420471589e-05;3.10043372564382e-05;3.24083420975188e-05;3.38425133444917e-05;3.53065457890846e-05;3.68005715186669e-05;3.83245389358488e-05;3.98782653882510e-05;4.14616952164644e-05;4.30748188780442e-05;4.47175650478249e-05;4.63898467993801e-05;4.80915319034599e-05;4.98226971677455e-05;5.15831512232700e-05;5.33728759952781e-05;5.51918321522111e-05;5.70398932535181e-05;5.89170846351397e-05;6.08231974387505e-05;6.27583343299031e-05;6.47222792528885e-05;6.67151027174107e-05;6.87366342110330e-05;7.07868142156144e-05;7.28656406545752e-05;7.49730151333333e-05;7.71088514536503e-05;7.92731023028393e-05;8.14657486020363e-05;8.36865503218298e-05;8.59356086301121e-05;8.82128718954112e-05;9.05181140965237e-05;9.28514356830927e-05;9.52127195198000e-05;9.76018480893283e-05;0.000100018833389875;0.000102463466950951;0.000104935888582002;0.000107435853506894;0.000109963427809404;0.000112518416931519;0.000115100878417533;0.000117710619245428;0.000120347656556815;0.000123011958360559;0.000125703392813205;0.000128421904837227;0.000131167359042717;0.000133939843777353;0.000136739170060611;0.000139565385326211;0.000142418283302892;0.000145297878316896;0.000148204149921329;0.000151136894680200;0.000154096147220371;];
 
  
  %initial conditions
  cVec = zeros(2*xPts,1);

  rVec = [-1.87224610995396e-16;-3.30786730057164e-12;-2.24724099531369e-11;4.01904154996178e-11;1.05515883490941e-10;9.30372762880140e-13;-2.22787146871473e-11;-9.01734922354510e-11;-1.10745807543975e-10;-7.46930791482446e-11;-7.04974264389784e-11;4.31467284250637e-11;1.12703422047960e-10;-4.36776668367891e-11;-9.80526426247680e-11;-1.69498941588624e-10;-1.50010671961555e-10;7.19193732481383e-11;4.90354184123562e-11;-3.73919643620498e-11;-2.30534115153733e-11;4.08012112865317e-11;1.97292880329322e-10;2.63783826352044e-10;-8.94713852160091e-11;1.32202516728991e-10;-6.02230789922050e-11;-2.15259724650231e-11;9.98554430331824e-11;6.13354586033715e-11;6.44512198060172e-11;-5.18601791827007e-11;-1.16298985686995e-10;7.65217226064187e-11;9.86204866143710e-11;-8.29445518345240e-12;-6.78931034114004e-11;3.02471279902024e-11;-2.08068015562172e-10;6.00770095563613e-11;4.54693481217322e-11;-4.80032759419061e-11;7.73088802176469e-11;-2.23744500574621e-11;-1.85586587330561e-11;2.01595816862882e-11;-1.67240320984169e-12;7.20655447209589e-12;-3.71625537417750e-11;-3.12902323409275e-11;-3.76327356863687e-11;-6.05724400758873e-12;2.23870513056321e-11;-1.81995921154712e-11;-1.49065678420671e-11;-1.48857884569101e-11;-8.20683636066494e-12;-1.20589816344358e-11;-1.68815471584399e-11;-5.10757517009630e-12;-2.19640629468831e-11;7.69523767624181e-13;-1.19783015971320e-11;-2.75028979625225e-12;8.88611580152077e-12;9.17275587083988e-13;7.01488476681500e-12;1.58665387686024e-11;1.49453431401436e-11;-6.52592368165228e-12;8.57727127170915e-12;5.94901797094927e-12;3.14754103857007e-12;-7.39040989863882e-13;-2.02692593612164e-12;-2.98649061210299e-12;-3.07789362573307e-12;-7.44722642928575e-13;-8.34349950240565e-13;1.36288094525205e-12;8.85102971817653e-13;2.06962919024889e-13;2.92512325125727e-13;-1.95851126115476e-12;5.33951572193048e-13;1.53728614039894e-12;4.03221510000829e-13;3.52591220109666e-14;1.81264291770900e-13;3.55048022232518e-14;1.17830658424856e-13;-6.71181860090186e-14;3.42939652370200e-14;-3.32487943427440e-14;1.09230116271397e-14;-2.66149949301742e-14;-2.48412401759879e-15;-6.27536217434610e-15;-3.33001855257198e-15;3.99680288865056e-15;0;-1.69909023257165e-18;-0.00180194456917740;-0.00359996877946229;-0.00539030161897191;-0.00716931304831092;-0.00893350624032737;-0.0106813283706023;-0.0124112815817571;-0.0141219201157762;-0.0158118476425007;-0.0174797147671106;-0.0191242167027280;-0.0207440910948473;-0.0223381159839777;-0.0239051078978705;-0.0254439200609780;-0.0269534407135857;-0.0284325915327688;-0.0298803261455125;-0.0312956287302472;-0.0326775126970926;-0.0340250194452503;-0.0353372171872712;-0.0366131998406438;-0.0378520859779137;-0.0390530178343277;-0.0402151603671504;-0.0413377003648335;-0.0424198456022618;-0.0434608240380305;-0.0444598830532369;-0.0454162887277894;-0.0463293251515870;-0.0471982937702556;-0.0480225127609642;-0.0488013164389084;-0.0495340546907948;-0.0502200924348154;-0.0508588091053545;-0.0514495981609144;-0.0519918666138635;-0.0524850345814093;-0.0529285348557132;-0.0533218124929829;-0.0536643244197531;-0.0539555390562795;-0.0541949359553346;-0.0543820054561420;-0.0545162483525224;-0.0545971755745197;-0.0546243078826509;-0.0545971755745534;-0.0545162483526045;-0.0543820054562457;-0.0541949359554697;-0.0539555390564727;-0.0536643244199783;-0.0533218124932349;-0.0529285348560089;-0.0524850345817141;-0.0519918666141957;-0.0514495981612692;-0.0508588091057694;-0.0502200924352542;-0.0495340546912524;-0.0488013164393973;-0.0480225127614573;-0.0471982937707453;-0.0463293251521281;-0.0454162887283069;-0.0444598830538092;-0.0434608240385718;-0.0424198456028313;-0.0413377003654400;-0.0402151603676662;-0.0390530178348595;-0.0378520859784781;-0.0366131998411603;-0.0353372171877931;-0.0340250194457342;-0.0326775126976656;-0.0312956287307328;-0.0298803261461042;-0.0284325915333013;-0.0269534407141305;-0.0254439200614265;-0.0239051078983431;-0.0223381159844076;-0.0207440910951981;-0.0191242167031094;-0.0174797147674447;-0.0158118476428555;-0.0141219201160822;-0.0124112815819589;-0.0106813283707926;-0.00893350624051167;-0.00716931304842916;-0.00539030161909891;-0.00359996877953503;-0.00180194456920927;0;];

  
  % Set force function.
  q    = zeros(xPts,1); % No force on most of it.
  fk = 10;
  q(5) = 1*fk;%*dx^4;
  q(end-4) = 1*fk;%*dx^4;
  q(mid) = -2*fk;%*dx^4;

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
            (1E-4/dx^2)*d2Coeffs2(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
    dOne = dOne + ...
            d1Coeffs1(count + 5) * diag( ones( 1, xPts - abs(count) ), count);
  end
  dTwo(1,:) = 0;
  dTwo(1,1) = 0;
  dTwo(end,:) = 0;
  dTwo(end,end) = 0;
  
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 = -100*dt/m*EBOp;                                                % map v to v stepping backwards
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
  while Y(5) <= 0.1;
     
    dw = cVec(xPts+2:end) - cVec(xPts+1:end-1) - deflection(2:end) + deflection(1:end-1);
    dphi = atan(dw/dx);
    
    for count2 = 1:mid - 1
      cphi(mid - count2) = cphi(mid - count2 + 1) + dphi(mid - count2);
      cphi(mid + count2) = cphi(mid + count2 - 1) + dphi(mid + count2 -1);
      
      X(mid - count2) = X(mid - count2 + 1) - dx*cos(cphi(mid - count2));
      X(mid + count2) = X(mid + count2 - 1) + dx*cos(cphi(mid + count2));
      
      Y(mid - count2) = Y(mid - count2 + 1) - dx*sin(cphi(mid - count2));
      Y(mid + count2) = Y(mid + count2 - 1) + dx*sin(cphi(mid + count2));   
    end
    
    if fap == 1
          woot = 0.1;
          sl = X(end - 4);
          sy0 = Y(5);
    end
    
    sdx = X(end - 4);        %string x component
    sa = acos(sdx/sl);              %angle of string
    sdy = sl*sin(sa);             %string y component

   F = [q/m*dt; zeros(xPts,1)]*sin(-cphi(5)+sa)/(2*sin(sa)+woot);
   F(mid) = q(mid)/m*dt;
   woot = 0;
   
    cVec = Tcrank * (cVec + F*fap/fin);                                             %CRANK
    cVec(xPts+1:2*xPts) = cVec(xPts+1:2*xPts) - cVec(mid+xPts);             %move back
    cVec(1:xPts) = cVec(1:xPts) - cVec(mid); 
    
        if mod(count, fin/10000) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        printF = F(5)
        sdx = sdx
        sa = sa
        hold off;
        figure(1)
        plot([Y(5)+sdy; Y(5); Y; Y(end-4); Y(5)+sdy],[0; X(5); X; X(end-4); 0]')
        hold on
        axis equal
        axis([-0.5 1 -1 1])
        end
    
    fap = fap + 1
  end
  cVec(1:xPts) = 0;
        figure(376)
        plot([Y(5)+sdy; Y(5); Y; Y(end-4); Y(5)+sdy],[0; X(5); X; X(end-4); 0]')
        hold on
        axis equal
        axis([-0.2 1 -1 1])
%CALCULATION LOOP----------------------------------------------------------
  q=q*1.1;
  hold off;
  for count = 1:fin;
     
    dw = cVec(xPts+2:end) - cVec(xPts+1:end-1) - deflection(2:end) + deflection(1:end-1);
    dphi = atan(dw/dx);
    
    for count2 = 1:mid - 1
      cphi(mid - count2) = cphi(mid - count2 + 1) + dphi(mid - count2);
      cphi(mid + count2) = cphi(mid + count2 - 1) + dphi(mid + count2 -1);
      
      X(mid - count2) = X(mid - count2 + 1) - dx*cos(cphi(mid - count2));
      X(mid + count2) = X(mid + count2 - 1) + dx*cos(cphi(mid + count2));
      
      Y(mid - count2) = Y(mid - count2 + 1) - dx*sin(cphi(mid - count2));
      Y(mid + count2) = Y(mid + count2 - 1) + dx*sin(cphi(mid + count2));   
    end
    
    if count == 1
          woot = 0.01;
          sl = X(end - 4);
          sy0 = Y(5);
    end
    
%     if count >=500
%         q = q*0;
%     end
    %q = q*(1+0.01*(sin(0.1*count)));
    
    sdx = X(end - 4);        %string x component
    sa = acos(sdx/sl);              %angle of string
    sdy = sl*sin(sa);             %string y component
    

   F = [q/m*dt; zeros(xPts,1)]*sin(-cphi(5)+sa)/(2*sin(sa)+woot);
   F(mid) = q(mid)/m*dt;
   woot = 0;
   
    cVec = Tcrank * (cVec + F);                                             %CRANK
    cVec(xPts+1:2*xPts) = cVec(xPts+1:2*xPts) - cVec(mid+xPts);             %move back
    cVec(1:xPts) = cVec(1:xPts) - cVec(mid);                                %decelerate
    
    %cVec = setBoundaries(cVec,xPts);
    
    if mod(count, fin/10000) == 0
        pause(0.01)
        clc
        pc = count*100/fin
        printF = F(5)
        sdx = sdx
        sa = sa
        hold off;
       figure(1)
%         plot(cVec(xPts+1:2*xPts));
%         axis([0 xPts -1E-3 1E-2])
%         figure(2)
%         plot(F(1:xPts))
%         figure(3)
%         plot(dphi)
%         axis([0 xPts -0.1 0.1])
%          figure(4)
%         plot(cphi)
%         axis([0 xPts -0.1 0.1])
        figure(1)
        plot([Y(5)+sdy; Y(5); Y; Y(end-4); Y(5)+sdy],[0; X(5); X; X(end-4); 0]')
        hold on
        axis equal
        axis([-0.5 1 -1 1])
    end
  end
