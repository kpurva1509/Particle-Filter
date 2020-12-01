function initialize_new_landmark(z, R)

global Param;
global State;



range = z(1);
bearing = z(2);
identity = z(3);

robot_x_bar = State.Ekf.mu_bar(1);
robot_y_bar = State.Ekf.mu_bar(2);
robot_theta_bar = State.Ekf.mu_bar(3);

map_x = robot_x_bar + range*cos(bearing + robot_theta_bar);
map_y = robot_y_bar + range*sin(minimizedAngle(bearing + robot_theta_bar));

variance = eye(2,2)*10000;

if length(State.Ekf.iM)==0
    map_x_index = State.Ekf.iR(end)+1;
else
    map_x_index = State.Ekf.iM(end)+1;
end
map_y_index = map_x_index + 1;

% Update iM matrix of state with new indices
indices_matrix = [map_x_index, map_y_index];
State.Ekf.iM = [State.Ekf.iM, indices_matrix];
State.Ekf.iL{end+1} = indices_matrix;
State.Ekf.sL(end+1) = identity;
State.Ekf.nL = State.Ekf.nL + 1;

% Update mu and sigma matrices
State.Ekf.mu(map_x_index) = map_x;
State.Ekf.mu(map_y_index) = map_y;
State.Ekf.mu_bar(map_x_index) = map_x;
State.Ekf.mu_bar(map_y_index) = map_y;
State.Ekf.Sigma = blkdiag(State.Ekf.Sigma, variance);
State.Ekf.sigma_bar = blkdiag(State.Ekf.sigma_bar, variance);

