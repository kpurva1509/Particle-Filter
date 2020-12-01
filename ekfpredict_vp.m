function ekfpredict_vp(u, delta_t)
% EKF-SLAM prediction for Victoria Park process model

global Param;
global State;

% Update mu using motion model
x_tm1 = State.Ekf.mu(1);
y_tm1 = State.Ekf.mu(2);
theta_tm1 = State.Ekf.mu(3);

a = Param.a;
b = Param.b;
L = Param.L;
H = Param.H;

ve = u(1);
alpha = u(2);
phi = theta_tm1;

vc = ve*1/(1 - tan(alpha)*H/L);

x_t = x_tm1 + delta_t*(vc*cos(phi) - vc/L*tan(alpha)*(a*sin(phi) + b*cos(phi)));
y_t = y_tm1 + delta_t*(vc*sin(phi) + vc/L*tan(alpha)*(a*cos(phi) - b*sin(phi)));
theta_t = minimizedAngle(phi + delta_t*vc/L*tan(alpha));

mu_bar = [x_t; y_t; theta_t];

N = State.Ekf.nL;

F_x = eye(3, 3+2*N);

G_t_3x3(:,1) = [1; 0; 0];
G_t_3x3(:,2) = [0; 1; 0];
G_t_3x3(:,3) = [delta_t*(-vc*sin(phi) - vc/L*tan(alpha)*(a*cos(phi) - b*sin(phi)));
                delta_t*(vc*cos(phi) + vc/L*tan(alpha)*(-a*sin(phi) - b*cos(phi)));
                1];

G_t = eye(3+2*N, 3+2*N) + F_x'*G_t_3x3*F_x;

V_t(:, 1) = [delta_t*(cos(phi) - tan(alpha)/L*(a*sin(phi) + b*cos(phi)));
             delta_t*(sin(phi) + tan(alpha)/L*(a*cos(phi) - b*sin(phi)));
             delta_t*tan(alpha)/L];
V_t(:, 2) = [delta_t*(-vc/L*sec(alpha)*sec(alpha)*(a*sin(phi) + b*cos(phi)));
             delta_t*(vc/L*sec(alpha)*sec(alpha)*(a*cos(phi) - b*sin(phi)));
             delta_t*vc/L*sec(alpha)*sec(alpha)];
   
M = Param.Qu;

R_t = V_t*M*V_t';

sigma_bar = G_t*State.Ekf.Sigma*G_t' + F_x'*R_t*F_x;

State.Ekf.mu_bar(1:3) = mu_bar(1:3);
State.Ekf.sigma_bar = sigma_bar;
