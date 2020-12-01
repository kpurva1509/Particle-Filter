function ekfpredict_sim(u)
% EKF-SLAM prediction for simulator process model

global Param;
global State;

% Update mu using motion model
prev_mu = State.Ekf.mu;
prev_sigma = State.Ekf.Sigma;
mu_bar = prediction(prev_mu, u);
mu_bar(3) = minimizedAngle(mu_bar(3));

% Update sigma using motion model
x = mu_bar(1);
y = mu_bar(2);
theta = mu_bar(3);

N = State.Ekf.nL;

F_x = eye(3, 3+2*N);

drot1 = minimizedAngle(u(1));
dtran = u(2);
drot2 = minimizedAngle(u(3));

G_t_3x3 = [0, 0, -dtran*sin(minimizedAngle(theta + drot1));
       0, 0,  dtran*cos(minimizedAngle(theta + drot1));
       0, 0,  0];

G_t = eye(3+2*N, 3+2*N) + F_x'*G_t_3x3*F_x;

V_t = [ -dtran*sin(minimizedAngle(theta + drot1)), cos(minimizedAngle(theta + drot1)), 0;
        dtran*cos(minimizedAngle(theta + drot1)), sin(minimizedAngle(theta + drot1)), 0;
        1,  0, 1];
   
M = diag([Param.alphas(1)*drot1^2+Param.alphas(2)*dtran^2, Param.alphas(3)*dtran^2+Param.alphas(4)*(drot1^2+drot2^2), Param.alphas(1)*drot2^2+Param.alphas(2)*dtran^2]);

R_t = V_t*M*V_t';

sigma_bar = G_t*prev_sigma*G_t' + F_x'*R_t*F_x;

State.Ekf.mu_bar = mu_bar;
State.Ekf.sigma_bar = sigma_bar;
