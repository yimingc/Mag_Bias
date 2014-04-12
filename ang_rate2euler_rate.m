% This function tranform angular rate omega_gb_b to euler angle rate
% i.e. the derivative of roll, pitch, yaw.
% Ref: Jay Farrell's Aided Navigation Chap 2.7
function euler_rate = ang_rate2euler_rate( ang_rate )
roll = ang_rate(1);
pitch = ang_rate(2);
yaw  = ang_rate(3);
Omega_E = [ 1  sin(roll)*tan(pitch)  cos(roll)*tan(pitch) 
               0        cos(roll)              -sin(pitch)          
               0  sin(roll)/cos(pitch)  cos(roll)/cos(pitch) ];
euler_rate = Omega_E*ang_rate;     