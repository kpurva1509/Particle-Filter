function varargout = runsim(stepsOrData, pauseLen)

global Param;
global Data;
global State;
global determinants;
determinants = [];

if ~exist('pauseLen','var')
    pauseLen = 0.3; % seconds
end

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
Param.initialStateMean = [180 50 0]';

% max number of landmark observations per timestep
Param.maxObs = 6;

% number of landmarks per sideline of field (minimum is 3)
Param.nLandmarksPerSide = 6;

% Motion noise (in odometry space, see p.134 in book).
Param.alphas = [0.05 0.001 0.05 0.01].^2; % std of noise proportional to alphas

% Standard deviation of Gaussian sensor noise (independent of distance)
Param.beta = [10, deg2rad(10)]; % [cm, rad]
Param.R = diag(Param.beta.^2);

% Step size between filter updates, can be less than 1.
Param.deltaT=0.1; % [s]
if isscalar(stepsOrData)
    % Generate a data set of motion and sensor info consistent with
    % noise models.
    numSteps = stepsOrData;
    Data = generateScript(Param.initialStateMean, numSteps, Param.maxObs, Param.alphas, Param.beta, Param.deltaT);
    save('current_run.mat', 'Data');
    
  % Generate a data set of motion and sensor info consistent with
    %load('current_run.mat')
else
    % use a user supplied data set from a previous run
    Data = stepsOrData;
    numSteps = size(Data, 1);
    global FIELDINFO;
    FIELDINFO = getfieldinfo;
end
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = Param.initialStateMean;
State.Ekf.Sigma = eye(3)*0.001;

for t = 1:numSteps
    plotsim(t);

    %=================================================
    % data available to your filter at this time step
    %=================================================
    u = getControl(t);
    z = getObservations(t);


    %=================================================
    %TODO: update your filter here based upon the
    %      motionCommand and observation
    %=================================================
    ekfpredict_sim(u);
    ekfupdate(z);
    State.Ekf.t = t;

    %=================================================
    %TODO: plot and evaluate filter results here
    %=================================================
    plotstate(z);

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

plot_corr(20)
plot_det(30)

if strcmp(lower(Param.dataAssociation), 'nn')
    figure;
    lengths_seq = 1:length(State.Ekf.nn_associations);
    scatter(lengths_seq, State.Ekf.nn_associations,'r*')
    hold on; 
    scatter(lengths_seq, State.Ekf.actual_associations,'b'),legend('Predicted','Actual');
end

if nargout >= 1
    varargout{1} = Data;
end
    
%==========================================================================
function u = getControl(t)
global Data;
% noisefree control command
u = Data.noisefreeControl(:,t);  % 3x1 [drot1; dtrans; drot2]


%==========================================================================
function z = getObservations(t)
global Data;
% noisy observations
z = Data.realObservation(:,:,t); % 3xn [range; bearing; landmark id]
ii = find(~isnan(z(1,:)));
z = z(:,ii);

%==========================================================================
function plotsim(t)
global Data;

%--------------------------------------------------------------
% Graphics
%--------------------------------------------------------------

NOISEFREE_PATH_COL = 'green';
ACTUAL_PATH_COL = 'blue';

NOISEFREE_BEARING_COLOR = 'cyan';
OBSERVED_BEARING_COLOR = 'red';

GLOBAL_FIGURE = 1;

%=================================================
% data *not* available to your filter, i.e., known
% only by the simulator, useful for making error plots
%=================================================
% actual position (i.e., ground truth)
x = Data.Sim.realRobot(1,t);
y = Data.Sim.realRobot(2,t);
theta = Data.Sim.realRobot(3,t);

% real observation
observation = Data.realObservation(:,:,t);

% noisefree observation
noisefreeObservation = Data.Sim.noisefreeObservation(:,:,t);

%=================================================
% graphics
%=================================================
figure(GLOBAL_FIGURE); clf; hold on; plotfield(observation(3,:));

% draw actual path (i.e., ground truth)//vicpark/viewLsr.m
plot(Data.Sim.realRobot(1,1:t), Data.Sim.realRobot(2,1:t), 'Color', ACTUAL_PATH_COL);
plotrobot( x, y, theta, 'black', 1, ACTUAL_PATH_COL);

% draw noise free motion command path
plot(Data.Sim.noisefreeRobot(1,1:t), Data.Sim.noisefreeRobot(2,1:t), 'Color', NOISEFREE_PATH_COL);
plot(Data.Sim.noisefreeRobot(1,t), Data.Sim.noisefreeRobot(2,t), '*', 'Color', NOISEFREE_PATH_COL);

for k=1:size(observation,2)
    rng = Data.Sim.noisefreeObservation(1,k,t);
    ang = Data.Sim.noisefreeObservation(2,k,t);
    noisy_rng = observation(1,k);
    noisy_ang = observation(2,k);

    % indicate observed range and angle relative to actual position
    plot([x x+cos(theta+noisy_ang)*noisy_rng], [y y+sin(theta+noisy_ang)*noisy_rng], 'Color', OBSERVED_BEARING_COLOR);

    % indicate ideal noise-free range and angle relative to actual position
    plot([x x+cos(theta+ang)*rng], [y y+sin(theta+ang)*rng], 'Color', NOISEFREE_BEARING_COLOR);
end

function plot_corr(fig)
    global State;
    x = corr(State.Ekf.Sigma);
    idx = [1,4:2:18];
    x_c = x(idx);
    x_c = x(idx, :);
    x_c = x_c(:,idx);
    figure;imagesc(x_c);title('correlation of x-coordinates of map with robot pose and x-coordinates');
    colorbar;
    
function plot_det(fig)
    global determinants
    figure; 
    for i = 1:8
        hold on; 
        plot(determinants(i,:));
    end
    legend('marker id 3','marker id 4', ...
    'marker id 5','marker id 6','marker id 7',...
    'marker id 8','marker id 1','marker id 2')

