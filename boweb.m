  clear
  clf
  xPts = 101;                                                               % Number of x points. Odd
  L = 1.5;
  dx = L/xPts;
  dt   = 0.00001;                                                              % Time step.
  m    = 1/1000;                                                             % Mass density.
  mid = (xPts+1)/2;
  
  %coupled PDE vector
  cVec = [zeros(2*xPts,1)];
  rVec = [-1.87224610995396e-16;-3.30786730057164e-12;-2.24724099531369e-11;4.01904154996178e-11;1.05515883490941e-10;9.30372762880140e-13;-2.22787146871473e-11;-9.01734922354510e-11;-1.10745807543975e-10;-7.46930791482446e-11;-7.04974264389784e-11;4.31467284250637e-11;1.12703422047960e-10;-4.36776668367891e-11;-9.80526426247680e-11;-1.69498941588624e-10;-1.50010671961555e-10;7.19193732481383e-11;4.90354184123562e-11;-3.73919643620498e-11;-2.30534115153733e-11;4.08012112865317e-11;1.97292880329322e-10;2.63783826352044e-10;-8.94713852160091e-11;1.32202516728991e-10;-6.02230789922050e-11;-2.15259724650231e-11;9.98554430331824e-11;6.13354586033715e-11;6.44512198060172e-11;-5.18601791827007e-11;-1.16298985686995e-10;7.65217226064187e-11;9.86204866143710e-11;-8.29445518345240e-12;-6.78931034114004e-11;3.02471279902024e-11;-2.08068015562172e-10;6.00770095563613e-11;4.54693481217322e-11;-4.80032759419061e-11;7.73088802176469e-11;-2.23744500574621e-11;-1.85586587330561e-11;2.01595816862882e-11;-1.67240320984169e-12;7.20655447209589e-12;-3.71625537417750e-11;-3.12902323409275e-11;-3.76327356863687e-11;-6.05724400758873e-12;2.23870513056321e-11;-1.81995921154712e-11;-1.49065678420671e-11;-1.48857884569101e-11;-8.20683636066494e-12;-1.20589816344358e-11;-1.68815471584399e-11;-5.10757517009630e-12;-2.19640629468831e-11;7.69523767624181e-13;-1.19783015971320e-11;-2.75028979625225e-12;8.88611580152077e-12;9.17275587083988e-13;7.01488476681500e-12;1.58665387686024e-11;1.49453431401436e-11;-6.52592368165228e-12;8.57727127170915e-12;5.94901797094927e-12;3.14754103857007e-12;-7.39040989863882e-13;-2.02692593612164e-12;-2.98649061210299e-12;-3.07789362573307e-12;-7.44722642928575e-13;-8.34349950240565e-13;1.36288094525205e-12;8.85102971817653e-13;2.06962919024889e-13;2.92512325125727e-13;-1.95851126115476e-12;5.33951572193048e-13;1.53728614039894e-12;4.03221510000829e-13;3.52591220109666e-14;1.81264291770900e-13;3.55048022232518e-14;1.17830658424856e-13;-6.71181860090186e-14;3.42939652370200e-14;-3.32487943427440e-14;1.09230116271397e-14;-2.66149949301742e-14;-2.48412401759879e-15;-6.27536217434610e-15;-3.33001855257198e-15;3.99680288865056e-15;0;-1.69909023257165e-18;-0.00180194456917740;-0.00359996877946229;-0.00539030161897191;-0.00716931304831092;-0.00893350624032737;-0.0106813283706023;-0.0124112815817571;-0.0141219201157762;-0.0158118476425007;-0.0174797147671106;-0.0191242167027280;-0.0207440910948473;-0.0223381159839777;-0.0239051078978705;-0.0254439200609780;-0.0269534407135857;-0.0284325915327688;-0.0298803261455125;-0.0312956287302472;-0.0326775126970926;-0.0340250194452503;-0.0353372171872712;-0.0366131998406438;-0.0378520859779137;-0.0390530178343277;-0.0402151603671504;-0.0413377003648335;-0.0424198456022618;-0.0434608240380305;-0.0444598830532369;-0.0454162887277894;-0.0463293251515870;-0.0471982937702556;-0.0480225127609642;-0.0488013164389084;-0.0495340546907948;-0.0502200924348154;-0.0508588091053545;-0.0514495981609144;-0.0519918666138635;-0.0524850345814093;-0.0529285348557132;-0.0533218124929829;-0.0536643244197531;-0.0539555390562795;-0.0541949359553346;-0.0543820054561420;-0.0545162483525224;-0.0545971755745197;-0.0546243078826509;-0.0545971755745534;-0.0545162483526045;-0.0543820054562457;-0.0541949359554697;-0.0539555390564727;-0.0536643244199783;-0.0533218124932349;-0.0529285348560089;-0.0524850345817141;-0.0519918666141957;-0.0514495981612692;-0.0508588091057694;-0.0502200924352542;-0.0495340546912524;-0.0488013164393973;-0.0480225127614573;-0.0471982937707453;-0.0463293251521281;-0.0454162887283069;-0.0444598830538092;-0.0434608240385718;-0.0424198456028313;-0.0413377003654400;-0.0402151603676662;-0.0390530178348595;-0.0378520859784781;-0.0366131998411603;-0.0353372171877931;-0.0340250194457342;-0.0326775126976656;-0.0312956287307328;-0.0298803261461042;-0.0284325915333013;-0.0269534407141305;-0.0254439200614265;-0.0239051078983431;-0.0223381159844076;-0.0207440910951981;-0.0191242167031094;-0.0174797147674447;-0.0158118476428555;-0.0141219201160822;-0.0124112815819589;-0.0106813283707926;-0.00893350624051167;-0.00716931304842916;-0.00539030161909891;-0.00359996877953503;-0.00180194456920927;0;];
   % Set force function.
  q    = zeros(xPts, 1);                                                   % No force on most of it.
  q(5) = 1000;
  q(end-4) = 1000;
  q((end+1)/2) = -2000;

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
  dTwo(1,:) = 0;
  dTwo(1,1) = -1;
  dTwo(end,:) = 0;
  dTwo(end,end) = -1;
  E = 1E9;
  I = [1:2/xPts:2];
  I = [I fliplr(I(1:end-1))];
  EI = diag(E.*I);
  
  %Euler-Bernoulli operator
  EBOp = dTwo*EI*dTwo;
  
  %Make HUGE MATRIX. Presently a huge operator such that 
  %Tfwd U(x,t) = U(x,t+1)
  M1 =  -10*dt/m*EBOp;                                      % map v to v stepping backwards
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
  fin = 10000;
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
    
    F = [q/m*dt; zeros(xPts,1)];%.*cos([cphi; zeros(xPts,1)]);
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
