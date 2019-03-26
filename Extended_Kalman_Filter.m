%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is written to filter out the measurement errors.
% The errors are filtered out using Extended Kalman Filter.
% The program generate states with synthetic noise and then sensors readings are used to estimate these states back.
% The sensor model is designed s.t. we have full state observability.
% - Prabhjeet Singh Arora
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Generating Synthetic Readings %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt   = .1;                      % Timestep for discrete sensor reading/output discrete time 0.1 sec
tf   = 240*60;                  % Total time 4hr - 240*60 sec
t    = 0:dt:tf;                 % Discrete time scale
m    = length(t);               % Length of the time scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% System of Manuevering Target %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A  =  diag([1 1 1 1],2);                            % Dynamic system model
B  = [zeros(4,2);eye(2)];                           % Control model
G  = [zeros(2);eye(2);zeros(2)];                    % Error model
Gn = [zeros(4,2);eye(2)];                           % Error model (for estimation for constant acceleration model)
%%%%% Continuous to Discrete (if discrete model is being used) %%%%%
[phi,gamw]  = c2d(A,G ,dt);                                                           % Error - dynamic model discretization
[phi,gamwn] = c2d(A,Gn,dt);                                                           % Error (acceleration) - dynamic model discretization
[phi,gam]   = c2d(A,B ,dt);                                                           % Control - dynamic model discretization

%%%%%%%%%%%%% Initial State %%%%%%%%%%%%%%%%%
x0 =  [10 15 11.11 8.33 0 0]';

%%%%%%%%%%%%% Control Input %%%%%%%%%%%%%%%%%
acc = 7.7160e-4;            % 10 km/hr^2 acceleration - 7.7160*10^-4 m/s^2
u   = zeros(m,2);
u((m-1)/4,1) =  (acc*dt);
u((m-1)/2,1) = -(acc*dt);
u((m-1)/4,2) =  (acc*dt);
u((m-1)/2,2) = -(acc*dt);

%%%%%%%%%%%%%%%%%% Generating truth model and error model and Estimation %%%%%%%%%%%%%%%%%
%%%%% Creating store variables to store generated synthetic states %%%%%
x_true = zeros(6,m);    % To store true states
x_true(:,1) = x0;
x_store = zeros(6,m);   % To store states with process noise
x_store(:,1) = x0;

%%%%% Error Covariances %%%%%
q = 10^(-2)*eye(2);     % Process noise error covariance
r = 10^(-2)*eye(2);     % Measurement noise error covariance

Qt = sqrt(q)*randn(2,m);
%%%%% Generating True States and Synthetic "Real" States with Process noise %%%%%
for i =2:m
    c = ode4(@fun,[t(i-1) t(i)],[x_true(:,i-1);u(i-1,:)';0;0],A,B,G);
    x_true(:,i) = c(end,1:6)';
    b = ode4(@fun,[t(i-1) t(i)],[x_store(:,i-1);u(i-1,:)';Qt(:,i-1)],A,B,G);
    x_store(:,i) = b(end,1:6)';
end

%%%%% Sensor Position at all time t %%%%%
%%%%%%%%%%%%%%%%%%%
% Change sensors to visualize the different types 3 sigma bounds.
% The way to determine position is through triangulation. When sensors and
% position of target are in line, the value is not determinable.
%%%%%%%%%%%%%%%%%%%
%%% Sensor 1 %%%
switch 2
    case 1
        X1_x = 10*(0.001*t)*1000;
        X1_y =  0*(0.005*t)*1000;
    case 2
        X1_x = -10*ones(1,m)*10^(4);
        X1_y =  zeros(1,m);
    case 3
        X1_x = -10*cos(0.001*t)*1000;
        X1_y =  5*sin(0.005*t)*1000;     
end
%%% Sensor 2 %%%
X2_x =  5*ones(1,m)*10^(4);
X2_y = -5*ones(1,m)*10^(4);

%%%%% Sensor True Reading for True States %%%%%
%%%%% Sensor 1 Reading for true state %%%%%
x1d_x = (X1_x - x_true(1,:)).^2;
x1d_y = (X1_y - x_true(2,:)).^2;
%%%%% Sensor 2 Reading for true state %%%%%
x2d_x = (X2_x - x_true(1,:)).^2;
x2d_y = (X2_y - x_true(2,:)).^2;
%%%%% Two sensors for triangulation %%%%%
%%%%% Sensor output %%%%%
y_true = [(x1d_x + x1d_y).^0.5;(x2d_x + x2d_y).^0.5];

%%%%% Sensor True Reading for Synthetic States %%%%%
%%%%% Sensor 1 Reading %%%%%
x1d_m_x = (X1_x - x_store(1,:)).^2;
x1d_m_y = (X1_y - x_store(2,:)).^2;
%%%%% Sensor 2 Reading %%%%%
x2d_m_x = (X2_x - x_store(1,:)).^2;
x2d_m_y = (X2_y - x_store(2,:)).^2;
%%%%% Two sensors for triangulation %%%%%
%%%%% Sensor output %%%%%
y_m = [(x1d_m_x+x1d_m_y).^0.5;(x2d_m_x+x2d_m_y).^0.5];
%%%%% Synthetic sensor output : Measurement noise implemented %%%%%
y_m = y_m + sqrt(r)*randn(2,m);

%%%%% Figure showing the synthetic position wrt true position %%%%%
figure(1)
measure = plot(x_store(1,:),x_store(2,:),'dr',x_true(1,:),x_true(2,:),'.-b',X1_x,X1_y,'.k',X2_x,X2_y,'ok','LineWidth',1.5);
%%%%% Plot properties %%%%%
%hold on
grid on
axis([-1 2 -1 2]*10^5)
xlabel('x-position (m)','FontSize',18)
ylabel('y-position (m)','FontSize',18)
title('Synthetic X,Y state and Sensor position','FontSize',20);
legend(measure,{'Synthetic "Real" Position','True Position','1st Sensor Position','2nd Sensor Position'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        Extended Kalman Filter        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial Conditions for Extended Kalman Filter %%%%%
xe0             = [10 15 11.11 8.33 0 0]';        % Initial State set for EKF this is also 'a priori'
xe              = xe0;
xe_store        = zeros(6,m);                     % Creating store variable to store estimated data
xe_store(:,1)   = xe0;
%yestore         = zeros(2,m);                    % If the output of sensor from estimated states needs to be stored

%%%%% Information Matrix Initialization %%%%%
pcov0           = 1*diag([1 1 1 1 1 1]);          % Creating Initial Information matrix / Covariance matrix
pcov            = pcov0;
pcov_store      = zeros(6,m);
pcov_store(:,1) = diag(pcov);                     % Creating the storage variable to store the diagonal of covariance matrix to get the 3 sigma bounds

%%%%% Error Covariance for Estimation %%%%%
Q = q;
R = r;

for ie = 2:m
    %%%%%%%%% PROPAGATION %%%%%%%%%%
    % Their are different options over the type of propagation change
    % switch number to access different options.
    switch 1
        %%% Discrete propagation %%%
        case 1
            xe   = phi*xe;
            switch 1
                case 1
                    pcov = phi*pcov*phi' + gamw*Q*gamw';
                case 2
                    pcov = phi*pcov*phi' + gamwn*Q*gamwn';
            end
        case 2
            xe   = phi*xe + gam*u(ie-1,:)';
            switch 1
                case 1
                    pcov = phi*pcov*phi' + gamw*Q*gamw';
                case 2
                    pcov = phi*pcov*phi' + gamwn*Q*gamwn';
            end    
        %%% Continuous propagation %%%
        case 3
            switch 1
                case 1
                    X    = ode4(@func,[t(ie-1) t(ie)],[xe;u(ie-1,:)';reshape(pcov,[],1)],A,B,G,Q);
                case 2
                    X    = ode4(@func,[t(ie-1) t(ie)],[xe;u(ie-1,:)';reshape(pcov,[],1)],A,B,Gn,Q);
            end
            xe   = X(end,1:6)';
            pcov = reshape(X(end,9:end)',6,6);
        case 4
            switch 2
                case 1
                    X    = ode4(@func,[t(ie-1) t(ie)],[xe;[0;0];reshape(pcov,[],1)],A,B,G,Q);
                case 2
                    X    = ode4(@func,[t(ie-1) t(ie)],[xe;[0;0];reshape(pcov,[],1)],A,B,Gn,Q);
            end
            xe   = X(end,1:6)';
            pcov = reshape(X(end,9:end)',6,6);
            
    end
    %%%%%%%%% GAIN %%%%%%%%%
    H = [(xe(1)-X1_x(ie))/((X1_x(ie)-xe(1))^2+(X1_y(ie)-xe(2))^2)^0.5, (xe(2)-X1_y(ie))/((X1_x(ie)-xe(1))^2+(X1_y(ie)-xe(2))^2)^0.5, 0, 0, 0, 0;
        (xe(1)-X2_x(ie))/((X2_x(ie)-xe(1))^2+(X2_y(ie)-xe(2))^2)^0.5, (xe(2)-X2_y(ie))/((X2_x(ie)-xe(1))^2+(X2_y(ie)-xe(2))^2)^0.5, 0, 0, 0, 0];
    K = pcov*H'/(H*pcov*H' + R);
    
    %%%%%%%%% UPDATE %%%%%%%%%
    ye   = [((X1_x(ie)-xe(1))^2+(X1_y(ie)-xe(2))^2)^0.5;((X2_x(ie)-xe(1))^2+(X2_y(ie)-xe(2))^2)^0.5];
    xe   = xe + K*(y_m(:,ie) - ye);
    pcov = (eye(6) - K*H)*pcov;
    
    %%%%%%%%% STORE %%%%%%%%%
    xe_store(:,ie)   = xe;
    pcov_store(:,ie) = diag(pcov);
end
sig3_kf = pcov_store.^(0.5)*3;      % To compute 3 sigma bounds

%%%%% Figure  %%%%%
figure(2)
subplot(2,1,1)
measuresub2 = plot(t,x_store(1,:),'dr',t,xe_store(1,:),'.b');
%%%%% Plot properties %%%%%
%hold on
grid on
%axis([0 tf -.5 .5])
xlabel('Time (sec)','FontSize',18);
ylabel('X-position (m)','FontSize',18);
title('Estimated states and Synthetic states','FontSize',20);
legend(measuresub2,{'Synthetic position','Estimated position'});
%%%%%%%%%%%%%%%%
subplot(2,1,2)
measuresub2 = plot(t,x_store(2,:),'dr',t,xe_store(2,:),'.b');
%%%%% Plot properties %%%%%
%hold on
grid on
%axis([0 tf -.5 .5])
xlabel('Time (sec)','FontSize',18);
ylabel('Y-position (m)','FontSize',18);
%title('Estimation Errors','FontSize',20);
legend(measuresub2,{'Synthetic position','Estimated position'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Figure  %%%%%
figure(3)
subplot(2,1,1)
measuresub3 = plot(t,xe_store(1,:)-x_store(1,:),'.b',t,sig3_kf(1,:),'r',t,-sig3_kf(1,:),'r');
%%%%% Plot properties %%%%%
%hold on
grid on
axis([0 tf -.5 .5])
ylabel('x-position error (m)','FontSize',18);
xlabel('Time (sec)','FontSize',18);
title('Estimation Errors','FontSize',20);
legend(measuresub3,{'Error','+3\sigma Bound','-3\sigma Bound'});
%%%%%%%%%%%%%%
subplot(2,1,2)
measuresub3 = plot(t,xe_store(2,:)-x_store(2,:),'.b',t,sig3_kf(2,:),'r',t,-sig3_kf(2,:),'r');
%%%%% Plot properties %%%%%
%hold on
grid on
axis([0 tf -.5 .5])
ylabel('y-position error (m)','FontSize',18);
xlabel('Time (sec)','FontSize',18);
%title('Estimation Errors','FontSize',20);
legend(measuresub3,{'Error','+3\sigma Bound','-3\sigma Bound'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% These functions are used for continuous propagation, using higher order
% ode solver.
function f1 = fun(~,x,A,B,G)
    xm = x(1:6);
    um = x(7:8);
    qm = x(9:10);
    xm_dot = A*xm+B*um+G*qm;
    f1 = [xm_dot;um;qm];
end
function f2 = func(~,x,A,B,G,Q)
    xa = x(1:6);
    ua= x(7:8);
    %xa_dot = A*xa;
    xa_dot = A*xa + B*ua;
    p = reshape(x(9:end),6,6);
    p_dot = A*p + p*A' + G*Q*G';
    f2 = [xa_dot;ua;reshape(p_dot,[],1)];
end