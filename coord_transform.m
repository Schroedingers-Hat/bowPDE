    dw1 = diff(cVec(xPts+1:end) + deflection);
    dw2 = diff(flipud(cVec(xPts+1:end) + deflection));
    phi1 = (atan(dw1/dx));
    phi2 = (atan(dw2/dx));
    phi = (phi1+phi2)/2;
    dphi = [0;diff(phi);0];
    somePhi = atan(dOne*(cVec(xPts+1:end) + deflection));
    for count2 = 1:mid - 1
      cphi(mid - count2) = cphi(mid - count2 + 1) - dphi(mid - count2);
      cphi(mid + count2) = cphi(mid + count2 - 1) + dphi(mid + count2);
      
      X(mid - count2) = X(mid - count2 + 1) - dx*cos(cphi(mid - count2));
      X(mid + count2) = X(mid + count2 - 1) + dx*cos(cphi(mid + count2));
      
      Y(mid - count2) = Y(mid - count2 + 1) - dx*sin(cphi(mid - count2));
      Y(mid + count2) = Y(mid + count2 - 1) + dx*sin(cphi(mid + count2));   
    end