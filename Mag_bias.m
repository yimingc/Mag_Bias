% This script is for mag sensor calibration.
% Refs: [1] Adaptive Estimation of Measurements Bias in Three-Dimensional
%           Field Sensors with Angular-Rate Sensors: Theory and Comparative
%           Experimental Evaluation.
% Notes: with conventional notation in navigation community
% 1. Measurement model of mag sensor
%    \tilde{x}(t) = x(t) + b + v,   
%    where x = R_g^b x_0 is the field vector in body frame
% 2. System dynamics, z = [x',b']
%    \dot{z}(t) = A(t)z(t) + w,  w~N(0,Q)
% 3. Measurement model
%    \tilde{y}(t) = Hz(t) + v,   v~N(0,R)
%    where H = [I3, I3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Constants      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_0 = [200, -40, 480]';                % ground truth of mag filed vector, in mG
bias   = [20, 120, 90]';               % ground truth of bias, in mG
ind_tm   = 1;                          % index for time in data struct
ind_mag  = 2:4;                        % index for mag meas in data struct
ind_gyro = 5:7;                        % index for gyro meas in data struct
ind_euler = 8:10;                      % index for euler angle in data struct
meas_size = length([ind_tm, ind_mag, ind_gyro]);
I3 = eye(3);
I6 = eye(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation options %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ite_num = 10;                          % MC iteration number
sim_time = 60;                         % simulation time length, in second
sigma_m = 1;                           % std dev of the mag meas noise, in mG
sigma_g = 5e-3;                        % std dev of the gyro meas noise, in rad/s
sim_freq = 100;                        % the freq to generate sim data, in Hz
ang_rate = pi/360;                       % magnitude of the generated angular rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation process %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate log memory
error_log = zeros(6, ite_num);
std_log   = zeros(6, ite_num);

for ind_sim = 1:ite_num
    % Generate sim data
    data_leng = sim_freq*sim_time;                               % number of data samples
    data = zeros(data_leng, meas_size);                          % allocate space for data
    dt = 1/sim_freq;                                             % dt
    data(:, ind_tm) = [dt:dt:sim_time]';                         % record time tags
    
    euler0 = rand(3,1)*pi/4;                                     % initial attitude in euler angle
    euler = euler0;
    omega = normrnd(0, ang_rate, data_leng, 3);                  % true angular rate
    omega_tilde = omega + normrnd(0, sigma_g, data_leng, 3);     % measured angular rate
    data(:, ind_gyro) = omega_tilde;                             % add to data struct
    for ind_data = 1:data_leng
        euler = euler + omega(ind_data,:)'*dt;
        data(ind_data, ind_euler) = limit_pi(euler);% add euler angle
        R_g2b = euler2R_g2b(limit_pi(euler));                              % get rotation matrix from global to body
        data(ind_data, ind_mag) = R_g2b*x_0 + bias;             % add mag meas
    end
    data(:, ind_mag) = data(:, ind_mag) + normrnd(0, sigma_m, data_leng, 3); % add noise to mag meas
    
    % Plot sim data, FIXME: try sphere plotting with mesh() function
    figure(1)
    plot(data(:, ind_tm),data(:, ind_euler));
    grid on;
    hold on;
    legend('roll', 'pitch', 'yaw');
    title('True Euler')
    xlabel('time (s)')
    ylabel('Angle (rad)')
    
    % Estimation
    z_hat_0 = [euler2R_g2b(data(1, ind_euler))*x_0+rand(3,1)*10;bias+rand(3,1)*10];% zeros(6,1);               
    z_hat = z_hat_0;                    % init estimate z_hat = [0,...,0]
    Q = blkdiag(I3*sigma_g,I3*sigma_g);                     
    R = I3*sigma_m;
    P_0 = I6*1;
    P_plus = P_0;
    H = [I3, I3];                       % measurement
    for ind_step = 1:data_leng
        omega_meas = data(ind_step, ind_gyro);
        y     = data(ind_step, ind_mag);
        % compute discrete system dynamics
        [Phi, Qd] = compute_system(omega_meas, Q, dt);
        % update
        z_hat_minus = Phi*z_hat;        
        P_minus = Phi*P_plus*Phi' + Qd;
        % innovation
        y_minus = H*z_hat_minus;
        dy = y' - y_minus;
        HPHtR = H*P_minus*H' + R;
        K = (P_minus*H')/(HPHtR);
        z_hat = z_hat_minus + K*dy;
        P_plus = (I6 - K*H)*P_minus;
    end
    error = z_hat(4:end) - bias;
    std   = diag(P_plus);
    disp(['Simulation ', num2str(ind_sim), ' done, the final estimation error is ', ...
        num2str(error')]);
    %pause();
end