function [R_n2b] = euler2R_g2b(x)
% convert Euler angles to a rotation matrix
% x = [roll, pitch, yaw]
c_r = cos(x(1));
s_r = sin(x(1));
c_p = cos(x(2));
s_p = sin(x(2));
c_y = cos(x(3));
s_y = sin(x(3));
R_n2b =[ c_y*c_p             s_y*c_p                 -s_p
    (-s_y*c_r+c_y*s_p*s_r) ( c_y*c_r+s_y*s_p*s_r)  (c_p*s_r)
    ( s_y*s_r+c_y*s_p*c_r) (-c_y*s_r+s_y*s_p*c_r)  (c_p*c_r)];  % eqn 2.43