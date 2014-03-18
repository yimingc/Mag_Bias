function [E] = R_b2g2euler(R)
% convert rotation matrix to Euler angles
E = zeros(3,1);
E(2) = atan2(-R(3,1),sqrt(R(1,1)^2+R(2,1)^2));
E(3) = atan2(R(2,1)/cos(E(2)),R(1,1)/cos(E(2)));
E(1) = atan2(R(3,2)/cos(E(2)),R(3,3)/cos(E(2)));