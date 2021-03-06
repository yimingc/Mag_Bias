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
%    where H = [I3, O3; hh, -hh] and 
%          hh = [z_hat(1)-z_hat(4),z_hat(2)-z_hat(5),
%          z_hat(3)-z_hat(6)]/norm_hat; which is the Jacobian
% 4. By toggling the flag 'with_norm_meas', you can choose whether include
%    the constant norm as a measurement.

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
norm_meas = norm(x_0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation options %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ite_num = 10;                          % MC iteration number
sim_time = 60;                         % simulation time length, in second
sigma_m = 1;                           % std dev of the mag meas noise, in mG
sigma_g = 5e-3;                        % std dev of the gyro meas noise, in rad/s
sim_freq = 100;                        % the freq to generate sim data, in Hz
ang_rate = pi/45;                      % magnitude of the generated angular rate
sigma_norm = 1e-0;                     % pick a std dev for the norm meas.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimation options %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_Q = 0.1;                           % Cov for system dynamics, see IV.B. [1]
var_x_init = 900;                      % variance for init x_0  estimae
var_bias_init = 900;                   % variance for init bias estimate
with_norm_meas = 1;                    % use |x-b| = 521.5 as a meas or not

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation process %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate log memory
data_leng = sim_freq*sim_time;                               % number of data samples
data = zeros(data_leng, meas_size);                          % allocate space for data
dt = 1/sim_freq;                                             % dt
error_log = zeros(data_leng+1, 3);                           % additional 1 for z_hat0
std_log   = zeros(data_leng+1, 3);                           % additional 1 for std0;

for ind_sim = 1:ite_num
    % Reset data struct
    data = data*0;
    error_log = error_log*0;
    std_log = std_log*0;
    
    % Generate sim data
    data(:, ind_tm) = [dt:dt:sim_time]';                         % record time tags
    euler0 = rand(3,1)*pi/4;                         % initial attitude in euler angle
    euler = euler0;
    R_g2b = euler2R_g2b(limit_pi(euler0));                       % get rotation matrix from global to body
    omega = normrnd(0, ang_rate, data_leng, 3);                  % true angular rate
    omega_tilde = omega + normrnd(0, sigma_g, data_leng, 3);     % measured angular rate
    data(:, ind_gyro) = omega_tilde;                     % add to data struct
    for ind_data = 1:data_leng
        ang_rate_temp = omega(ind_data,:)';                      
        euler_rate = ang_rate2euler_rate(ang_rate_temp);
        euler = euler + euler_rate*dt;
        data(ind_data, ind_euler) = limit_pi(euler);             % add euler angle
        R_g2b = (I3 - vcross(ang_rate_temp)*dt)*R_g2b;
        data(ind_data, ind_mag) = R_g2b*x_0 + bias;              % add mag meas
    end
    data(:, ind_mag) = data(:, ind_mag) + normrnd(0, sigma_m, data_leng, 3); % add noise to mag meas
    
    % Plot sim data, FIXME: try sphere plotting with mesh() function
    figure(1)
    plot(data(:, ind_tm), data(:, ind_euler));
    grid on;
    hold on;
    legend('roll', 'pitch', 'yaw');
    title('True Euler')
    xlabel('time (s)')
    ylabel('Angle (rad)')
    
    % Estimation
    z_hat_0 = [euler2R_g2b(euler0)*x_0+bias+rand(3,1)*sqrt(var_x_init);bias+rand(3,1)*sqrt(var_bias_init)];% zeros(6,1);
    %z_hat_0 = [euler2R_g2b(euler0)*x_0+bias;bias];
    z_hat = z_hat_0;                    % init estimate z_hat
    Q = var_Q*eye(6);                   % instead of Q = blkdiag(I3*sigma_g,I3*sigma_g),
    if with_norm_meas
        R = blkdiag(I3*sigma_m^2, sigma_norm^2);
    else
        R = I3*sigma_m^2;
    end
    P_0 = blkdiag(I3*var_x_init, I3*var_bias_init);
    P_plus = P_0;
    error_log(1,:) = z_hat(4:end) - bias;
    std = sqrt(diag(P_plus));
    std_log(1, :) = std(4:end);
    dx_minus = zeros(6,1);
    dx_plus = zeros(6,1);
    for ind_step = 1:data_leng
        omega_meas = data(ind_step, ind_gyro);
        y     = data(ind_step, ind_mag);
        if with_norm_meas
            norm_hat = norm(z_hat(1:3)-z_hat(4:6));        % estimated norm
            hh = [z_hat(1)-z_hat(4),z_hat(2)-z_hat(5), z_hat(3)-z_hat(6)]/norm_hat; % for the norm meas
            H = [I3, zeros(3,3)
                 hh, -hh];       % measurement model
        else
            H = [I3, zeros(3,3)];
        end
        % compute discrete system dynamics
        [Phi, Qd] = compute_system(omega_meas, Q, dt);
        Ob = obsv(Phi,H);
        unob = length(Phi)-rank(Ob);    % check observability
        % update
        z_hat_minus = Phi*z_hat;
        P_minus = Phi*P_plus*Phi' + Qd;
        % innovation
        y_minus = H*z_hat_minus;
        if with_norm_meas
            dy = [y';norm_meas] - y_minus;
        else
            dy = y' - y_minus;
        end
        HPHtR = H*P_minus*H' + R;
        K = (P_minus*H')/(HPHtR);
        z_hat = z_hat_minus + K*dy;
        P_plus = (I6 - K*H)*P_minus;
        error_log(ind_step+1, :) =  z_hat(4:end) - bias;
        std = sqrt(diag(P_plus));
        std_log(ind_step+1, :) = std(4:end);
    end
    % Result plotting and logging
    figure(2);
    subplot(3,1,1);
    title('Bias estimation error, with 3 sigma envelope')
    plot([0;data(:, ind_tm)], error_log(:,1), '*-');
    hold on; grid on; legend('Mag_x Est error');
    plot([0;data(:, ind_tm)], 3*std_log(:,1), 'r-.', [0;data(:, ind_tm)], -3*std_log(:,1), 'r-.');
    subplot(3,1,2);
    plot([0;data(:, ind_tm)], error_log(:,2), '*-');
    hold on; grid on; legend('Mag_y Est error');
    plot([0;data(:, ind_tm)], 3*std_log(:,2), 'r-.', [0;data(:, ind_tm)], -3*std_log(:,2), 'r-.');
    subplot(3,1,3);
    plot([0;data(:, ind_tm)], error_log(:,3), '*-');
    hold on; grid on; legend('Mag_z Est error');
    plot([0;data(:, ind_tm)], 3*std_log(:,3), 'r-.', [0;data(:, ind_tm)], -3*std_log(:,3), 'r-.');
    
    error = z_hat(4:end) - bias;
    norm(z_hat(1:3)-z_hat(4:end))
    std   = diag(P_plus);
    disp(['Simulation ', num2str(ind_sim), ' done, the final estimation error is ', ...
        num2str(error')]);
    break;
    %pause();
end