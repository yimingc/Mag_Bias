% This script is for mag sensor calibration.

% Notes: with conventional notation in navigation community
% 1. Measurement model of mag sensor
%    \tilde{x}(t) = x(t) + b + v,   
%    where x = R_g^b * x_0 is the field vector in body frame
% 2. System dynamics, z = [x',b']
%    \dot{z}(t) = A(t)z(t) + w,  w~N(0,Q)
% 3. Measurement model
%    \tilde{y}(t) = Hz(t) + v,   v~N(0,R)
%    where H = [I3, I3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Constants      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = bias_estimation()
clear
close all;

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
ite_num = 1;                                % MC iteration number
sim_time = 15;                            % simulation time length, in second
sigma_m = 1;                                % std dev of the mag meas noise, in mG
sigma_g = 5e-3;                             % std dev of the gyro meas noise, in rad/s
sim_freq = 100;                             % the freq to generate sim data, in Hz
ang_rate = pi/4;                           % magnitude of the generated angular rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation process %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate log memory
error_log = zeros(6, ite_num);
std_log   = zeros(6, ite_num);

% fname_DCM='.\INVN02191646PM.xlsx';
% Logdata=xlsread(fname_DCM); 
% gyroData(:,1) = Logdata(:,2)*8/131;
% gyroData(:,2) = Logdata(:,3)*8/131;
% gyroData(:,3) = Logdata(:,4)*8/131;
% omega = gyroData;

for ind_sim = 1:ite_num
    
    % Generate sim data
    data_leng = sim_freq*sim_time;                               % number of data samples
    data = zeros(data_leng, meas_size);                          % allocate space for data
    dt = 1/sim_freq;                                             % dt
    data(:, ind_tm) = [dt:dt:sim_time]';                         % record time tags
    
    %generating mag measurement and noise
    euler0 = rand(3,1)*pi/4;                                     % initial attitude in euler angle
    euler = euler0;
    omega = normrnd(0, ang_rate, data_leng, 3);                  % true angular rate
    
    for ind_data = 1:data_leng
        euler = euler + (omega(ind_data,:)*dt)';
        data(ind_data, ind_euler) = limit_pi(euler);             % add euler angle
        R_g2b = euler2R_g2b(limit_pi(euler));                    % get rotation matrix from global to body
        data(ind_data, ind_mag) = R_g2b*x_0 + bias;              % add mag meas
    end
    data(:, ind_mag) = data(:, ind_mag);% + normrnd(0, sigma_m, data_leng, 3); % add noise to mag meas
    
    %generating gyro measurement and noise
    omega_tilde = omega;% + normrnd(0, sigma_g, data_leng, 3);     % measured angular rate
    data(:, ind_gyro) = omega_tilde;                             % add to data struct
    
    % Plot sim data, FIXME: try sphere plotting with mesh() function
    figure(1)
    plot(data(:, ind_tm),data(:, ind_euler));
    grid on;
    hold on;
    legend('roll', 'pitch', 'yaw');
    title('True Euler')
    xlabel('time (s)')
    ylabel('Angle (rad)')
    
    %Plot mag data on 3D
    figure(2)
    [x,y,z] = sphere(30);colormap([1,1,1])
    surf(522*x+20, 522*y+120, 522*z+90) % where (a,b,c) is center of the sphere
    hold on;
    plot3(data(:, 2),data(:, 3),data(:, 4),'r.');
    axis equal;
    
    % Estimation
    %z_hat_0 = [euler2R_g2b(data(1, ind_euler))*x_0+rand(3,1)*10; rand(3,1)*10];% zeros(6,1);               
    
    z_hat_0 =  zeros(6,1);               
    z_hat = z_hat_0;                            % init estimate z_hat = [0,...,0]
    Q = blkdiag(I3*sigma_g,I3*sigma_g);                     
    R = I3*sigma_m;
    P_0 = I6*1;
    P_plus = P_0;
    H = [I3, zeros(3,3)];                       % measurement
    %H = [I3, I3];                       % measurement
    
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
        std(1:6,ind_step)   = diag(P_plus);
        
        z_hat_Store(1:3,ind_step)=z_hat(4:6);
        error(1:3,ind_step) = z_hat(4:end) - bias;        

    end

    [figure3] = figure(3);
    plot(z_hat_Store(1,:),'r.' );hold on
    plot(z_hat_Store(2,:),'g.' );hold on
    plot(z_hat_Store(3,:),'b.' );hold on
    addFigureTitle( figure3, 'Bias' );

    [figure4]= figure(4)
    subplot(3,1,1);
    plot(error(1,:),'r.');
    subplot(3,1,2);
    plot(error(2,:),'g.');
    subplot(3,1,3);
    plot(error(3,:),'b.');
    hold on
    addFigureTitle( figure4, 'Error of Bias' );
    
    [figure5]=figure(5)
    subplot(3,1,1);
    plot(std(1,:),'r.');
    subplot(3,1,2);
    plot(std(2,:),'g.');
    subplot(3,1,3);
    plot(std(3,:),'b.');
    hold on;
    addFigureTitle( figure5, 'Std Error of Bias' );

    [figure6] = figure(6);
    plot(data(:,2),'r.' );hold on
    plot(data(:,3),'g.' );hold on
    plot(data(:,4),'b.' );hold on; grid on
    addFigureTitle( figure6, 'mag' );

    [figure7] = figure(7);
    plot(data(:,5)*180/pi,'r.' );hold on
    plot(data(:,6)*180/pi,'g.' );hold on
    plot(data(:,7)*180/pi,'b.' );hold on;  grid on
    addFigureTitle( figure7, 'gyro' );

    [figure8] = figure(8);
    plot(data(:,8),'r.' );hold on
    plot(data(:,9),'g.' );hold on
    plot(data(:,10),'b.' );hold on;  grid on
    addFigureTitle( figure8, 'euler' );

%    error = z_hat(4:end) - bias;
%    std   = diag(P_plus);
    disp(['Simulation ', num2str(ind_sim), ' done, the final estimation error is ', ...
        num2str(error(1:3,ind_step)')]);
    %pause();
end

end

function [ output_args ] = addFigureTitle( figure1, str )
% Add a title near the top of a figure using the annotation function.
% This function is useful for labeling figures with multiple subplots.

annotation(figure1,'textbox',...
    [0.5 0.95 .1 0.05],...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','none',...
    'FitBoxToText','off',...
    'HorizontalAlignment', 'center', ...
    'Interpreter', 'none', ...
    'String',{str}, 'FitBoxToText','on');
end
