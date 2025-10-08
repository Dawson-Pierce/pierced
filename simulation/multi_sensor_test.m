%% Multi-Sensor Point Cloud simulations
%
% Author: Dawson Pierce

clear; clc; close all

addpath('../../brew/')

%% Import STL for TR obj

TR = stlread("drone_shape.stl");

%% Setup motion model / ETT PHD

target_motion = BREW.dynamics.Integrator_3D();

% generate measurements 

dt = 0.2;
tf = 10; 

t = 0:dt:tf;

measurements = {};

% initialize filter(s)
H = [eye(3) zeros(3, 3)];
R = 0.2 * eye(3); % Measurement noise
process_noise = 0.01 * eye(6); 
inner_filter = BREW.filters.TrajectoryGGIWEKF('dyn_obj',target_motion, ...
    'H',H,'process_noise',process_noise,'measurement_noise', R, 'L',20);

alpha = {10};
beta = {1}; 
mean = {[15; 15; 10; 0; 0; 0]}; 
covariance = {diag([10,10,10,2,2,2])}; 
IWdof = {10}; 
IWshape = {[1 0 0; 0 1 0; 0 0 1]}; 
weight = [1];

birth_model = BREW.distributions.TrajectoryGGIWMixture( ...
    'idx', {1}, ...
    'alphas', alpha, ...
    'betas', beta, ...
    'means', mean, ...
    'covariances', covariance, ...
    'IWdofs', IWdof, ...
    'IWshapes', IWshape, ...
    'weights', weight);

% phd = BREW.multi_target.PHD('filter',inner_filter, 'birth_model', birth_model,...
%     'prob_detection', 0.8, 'prob_survive', 0.9, 'max_terms',50, ...
%     'cluster_obj',BREW.clustering.DBSCAN_obj(3,5),'extract_threshold',0.5);

%% Setup sensor_station objs

sensors{1} = sensor_station('location',[0 20 0],'view_vector',[1 0 0],'grid_slope',0.03);
sensors{2} = sensor_station('location',[20 0 0],'view_vector',[0 1 0],'grid_slope',0.03);

% Setting up internal trackers this way so we can create copies since the
% classes for tracking inherit the handle class
for k = 1:length(sensors)
    sensors{k}.internal_tracker = BREW.multi_target.PHD( ...
        'filter',copy(inner_filter), ...
        'birth_model', copy(birth_model),...
    'prob_detection', 0.8, ...
    'prob_survive', 0.9, ...
    'max_terms',50, ...
    'cluster_obj',BREW.clustering.DBSCAN_obj(3,5), ...
    'extract_threshold',0.5);
end

T = [eye(3), [20; 20; 10]; zeros(1,3) 1];

points = cell(length(sensors),1); 

%% Simulation

x_true = [20; 20; 10; -1; -1; -0.1];

f = figure; 
subplot(1,2,1)
ax = axes;
axis equal; 
view(3); hold on

subplot(1,2,2)
ax2 = axes;
axis equal; 
view(3); hold on

est_mix = cell(length(sensors),1);
est_colors = {'m','c'};

gifFilename = 'PointCloudDroneTrackingPOC.gif';

for k = 1:length(t) 
    x_true = target_motion.propagateState(dt,x_true);
    T(1:3,4) = x_true(1:3);
    subplot(1,2,1)
    cla 
    for kk = 1:length(sensors)
        points{kk} = sensors{kk}.get_points(TR,T);
        scatter3(points{kk}(:,1),points{kk}(:,2),points{kk}(:,3),6,'*'); hold on
        sensors{kk}.plot()
        
        sensors{kk}.internal_tracker.predict(dt,{});
        sensors{kk}.internal_tracker.correct(dt,points{kk}');
        
        est_mix{kk} = sensors{kk}.internal_tracker.cleanup();
        
        est_mix{kk}.plot_distributions([1,2,3],'Color',est_colors{kk},'window_color','r','window_width',1.6);
    end

    xlim([0 30])
    ylim([0 30])
    zlim([-6 12])

    subplot(1,2,2)
    cla 
    for kk = 1:length(sensors)
        points{kk} = sensors{kk}.get_points(TR,T);
        scatter3(points{kk}(:,1),points{kk}(:,2),points{kk}(:,3),6,'*'); hold on
        sensors{kk}.plot()
    end

    xlim([0 30])
    ylim([0 30])
    zlim([-6 12])
    
    drawnow;
    
    frame = getframe(f);
    img   = frame2im(frame);
    
    % convert to an indexed image with a fixed 256‐color map
    [A, map] = rgb2ind(img, 256);
    
    % write to GIF: first frame creates the file, subsequent frames append
    if k == 1
        imwrite(A, map, gifFilename, 'gif', ...
                'LoopCount', Inf, ...      % make it loop forever
                'DelayTime', 0.1);         % seconds between frames
    else
        imwrite(A, map, gifFilename, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', 0.1);
    end

end

%%

figure;
for kk = 1:length(sensors) 
    scatter3(points{kk}(:,1),points{kk}(:,2),points{kk}(:,3),6,'*'); hold on
    sensors{kk}.plot()
end