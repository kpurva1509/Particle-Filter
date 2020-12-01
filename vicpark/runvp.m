function runvp(nSteps,pauseLen)
close all;

global Param;
global State;
global Data;

if ~exist('nSteps','var') || isempty(nSteps)
    nSteps = inf;
end

if ~exist('pauseLen','var')
    pauseLen = 0; % seconds
end

Data = load_vp_si();

makeVideo = 1;

if makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end

% Initalize Params
%===================================================
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = Param.Qf;
State.Ekf.mu_orig = State.Ekf.mu;

State.Ekf.mu_bar = State.Ekf.mu;
State.Ekf.sigma_bar = State.Ekf.Sigma;


global AAr;
AAr = [0:360]*pi/360;


figure(1); clf;
axis equal;

ci = 1; % control index
t = min(Data.Laser.time(1), Data.Control.time(1));

global Statistics;
Statistics = struct();
Statistics.delta_time_predict = [];
Statistics.delta_time_update = [];
Statistics.num_landmarks = [];


for k=1:min(nSteps, length(Data.Laser.time))
    
    tic;
    while (Data.Control.time(ci) < Data.Laser.time(k))
       % control available
       dt = Data.Control.time(ci) - t;
       t = Data.Control.time(ci);
       u = [Data.Control.ve(ci), Data.Control.alpha(ci)]';
       display(u);
       display(dt);
       display(k);
       ekfpredict_vp(u, dt);
       ci = ci+1;
    end
    delta_time = toc;
    Statistics.delta_time_predict(end+1) = delta_time;
    
    % observation available
    tic;
    dt = Data.Laser.time(k) - t;
    t = Data.Laser.time(k);
    z = detectTreesI16(Data.Laser.ranges(k,:));
    z_temp = z;
    z_temp(2, :) = z_temp(2, :) - pi/2;
    ekfupdate(z_temp);
    
    delta_time = toc;
    Statistics.delta_time_update(end+1) = delta_time;
    
    Statistics.num_landmarks(end+1) = State.Ekf.nL;
    
    doGraphics(z);
    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end
    
    if makeVideo
        F = getframe(gcf);
        switch votype
            case 'avifile'
                vo = addframe(vo, F);
            case 'VideoWriter'
                writeVideo(vo, F);
            otherwise
                error('unrecognized votype');
        end
    end    
end

display(State.Ekf)

% Plot the required statistics
num_steps = length(Statistics.num_landmarks);

figure;
plot(1:num_steps, Statistics.num_landmarks, '*');
xlabel('time-steps')
ylabel('number of detected landmarks')
ylim([5, 100]);

figure;
plot(1:num_steps, Statistics.delta_time_predict, '*');
xlabel('time-steps')
ylabel('predict CPU time consumed')

figure;
plot(1:num_steps, Statistics.delta_time_update, '*');
xlabel('time-steps')
ylabel('update CPU time consumed')

%==========================================================================
function doGraphics(z)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!

global Param;
global State;

% plot the robot and 3-sigma covariance ellipsoid
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on;

plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma, 'blue', 0, 'blue', 0, 3);

% restrict view to a bounding box around the current pose
BB=100;
axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);

mu = State.Ekf.mu;
Sigma = State.Ekf.Sigma;
robot_pos = mu(1:3);
robot_variance = Sigma(1:2, 1:2);

plotcov2d(mu(1), mu(2), Sigma(1:2, 1:2), 'b', 0, 0, 0, 3);

% project raw sensor detections in global frame using estimate pose
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr - pi/2);
    yl = yr + r*sin(b+tr - pi/2);
    plot([xr; xl], [yr; yl],'r',xl,yl,'r*');
end

for index = 1:State.Ekf.nL
    indices = State.Ekf.iL{index};
    map_x = mu(indices(1));
    map_y = mu(indices(2));
    map_sigma = Sigma(indices(1):indices(2), indices(1):indices(2));
    plotcov2d(map_x, map_y, map_sigma, 'blue', 0, 0, 0, 3);
end

hold off;

