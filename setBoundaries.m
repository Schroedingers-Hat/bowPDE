  function outVec = setBoundaries(Vec,xPts)

Vec((xPts+1)/2) = 0;
Vec((3*xPts+1)/2) = 0;

bw = 4;

for k = 1:bw
    Vec(k) = Vec(bw+1);
    Vec(xPts + 1 - k) = Vec(xPts - bw);
    
    Vec(xPts + bw + 1 - k) = Vec(xPts + bw + 2 - k) + Vec(xPts + bw + 1) - Vec(xPts + bw + 2);
    Vec(end - bw + k) = Vec(end - bw - 1 + k) + Vec(end - bw) - Vec(end - bw - 1);
end


outVec = Vec;