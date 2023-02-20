function onm_i = omeganm_i(r,dr,qnm_i_1,qnm_i_2,qnm_i_3)
% omega^{(nm)}_i
%
% Hai 11/27/22

dx = dr(1,:); dy = dr(2,:); dz = dr(3,:);
onm_i =  (dx.*r(3,:))'.*qnm_i_2-(dx.*r(2,:))'.*qnm_i_3 ...
        +(dy.*r(1,:))'.*qnm_i_3-(dy.*r(3,:))'.*qnm_i_1 ...
        +(dz.*r(2,:))'.*qnm_i_1-(dz.*r(1,:))'.*qnm_i_2;

end