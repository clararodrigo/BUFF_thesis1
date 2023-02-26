%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/Verasonics_GPU_Beamformer_CheeHau/src/gpuDAS'));
%BUFF
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/01 Save Fast'));
addpath(genpath('../buff/src'));
%FieldII
addpath(genpath("../field_ii"));
%save directories
addpath(genpath('RF'));
addpath(genpath('BF'));

field_init(0);
set_sampling(GlobalConfig().fs);


%% Setup
% Main volume
space = Box( ...
        [0, 0, 70*mm], ...              % Center
        [20*mm, 20*mm, 20*mm], ...      % Size
        [0, 0, 0] ...                   % Rotation
);


transducer = Matrix1024();
transducer.initial_load_apodizations('saved500CompRndApod.mat');
transducer.tx_aperture.excitation.n_cycles = 1;
transducer.tx_aperture.apply_waveforms();
transducer.tx_aperture.apply_delays();
transducer.tx_aperture.apply_apodization();
transducer.set_MI(0.05, space.center);

% apply steering delays (because its constant)
lambda = transducer.c/(transducer.f0);
delays = transducer.tx_aperture.calc_delays_diverging_pos([0,0,-32*lambda]);
for ee = 1:transducer.tx_aperture.n_elements
    transducer.tx_aperture.elements(ee).delay = delays(ee);
end
transducer.tx_aperture.apply_delays();
%% PLOT VIRTUAL SOURCES

figure;
ha = axes(); xlabel('X'); ylabel('Y'); zlabel('Z'); hold(ha, 'on');
xlim(ha, [-30e-3, 30e-3]);
ylim(ha, [-30e-3, 30e-3]);
zlim(ha, [-10e-3, 30e-3]);
grid on;
set(gca, 'ZDir','reverse')
transducer.tx_aperture.plot_aperture(ha);
% space.plot_skeleton(ha);
scatter3(ha, 0,0,-32*lambda,30,'black','filled')
scatter3(ha, 0,0,10*32*lambda,30,'black','filled')

delays_x = repmat([1:32],32,1);
scatter3(ha,delays_x,delays_x,delays)


%%
tube1 = PhantomTube(...
        [0, -0.2*mm, 70*mm], ...            % Center
        [200*um, 200*um, 30*mm], ...    % Size
        [0, 50, 90] ...                % Rotation
);
tube2 = PhantomTube(...
        [0, 0.2*mm, 70*mm], ...            % Center
        [200*um, 200*um, 30*mm], ...    % Size
        [0, -50, 90] ...                 % Rotation
);

% place bubbels
fps = 500;
dt = 1/fps;
dt = 0.01;
figure;
ha = axes(); xlabel('X'); ylabel('Y'); zlabel('Z'); hold(ha, 'on');
tube1.plot_skeleton(ha);
tube2.plot_skeleton(ha);

transducer.tx_aperture.plot_aperture(ha);
space.plot_skeleton(ha);

% pre-sim
tube1.update(1);
tube2.update(1);
view(3);
axis image

tube1.bub_scat.bub(6:end,:) = [];
tube1.bub_scat.pos(6:end,:) = [];
tube2.bub_scat.bub(6:end,:) = [];
tube2.bub_scat.pos(6:end,:) = [];

bub1_p = tube1.plot_bub_scat(ha);
bub2_p = tube2.plot_bub_scat(ha);

% pre-run
for fn = 1:UserSet.totalFrame
    tube1.update(dt);
    tube2.update(dt);
end


% Simulate
for f = 1:UserSet.totalFrame
    
%     waitbar(f/UserSet.totalFrame, ['Frame ',num2str(f),'/',num2str(UserSet.totalFrame)]);
    tube1.update(dt);
    tube2.update(dt);

    tube1.update_bub_scat(bub1_p);
    tube2.update_bub_scat(bub2_p);

    tot_bub(:,f) = tube1.bub_scat + tube2.bub_scat;
        
    fprintf('\n');
end

%% Pad signals 
max_t = max(arrayfun(@(x) max(bub_rf(x).time_vector), 1:numel(bub_rf)));

for ii = 1:numel(bub_rf)
	bub_rf(ii) = RFSignals(zeros(1024,2), 0) + bub_rf(ii) + RFSignals(zeros(1024,2), max_t);
end

%% Cat RF data for all apodizations
for f = 1:UserSet.totalFrame
    for ap = 1:UserSet.ap
        rf_data(:,:,f,ap) = bub_rf(f,ap).signals;
    end
end
rf_data = permute(rf_data, [2,1,3,4]);

% compounding
rf_data = sum(rf_data,4);

%% to beamform in gpu
filename  = [filename, '_dt', num2str(dt)];
save(['../RF/rf_', filename], 'rf_data', '-v7.3');
save(['../RF/rf_', filename,'bubbles'], 'tot_bub','-v7.3');
%%
% perform SVD
% svd_data = channel_svd(rf_data, 10, 50);

%% Beamforming
for i = 1:40:161
    % by 10-frame chunks
    [tmp, pixelMap] = beamform_sim_ple(rf_data(:,:,i:i+39), transducer, [0,0], -32*lambda);
    tmp = sum(tmp,5);
    tmp = permute(tmp, [3,2,1,4]);
    saveIQ_simple(tmp, ['../BF/bf_', filename,'_', num2str(i),'_',num2str(i+39)], pixelMap, UserSet)
end
%%


% function plot_res(img, pixelMap)
%     t = permute(temp_out, [3,2,1,4,5]);
%     t = abs(t(:,:,:,2));
%     p_x = squeeze(sum(t,1));
%     p_y = squeeze(sum(t,2));
%     p_z = squeeze(sum(t,3));
% 
%     figure;
%     subplot(131); imagesc(pixelMap.pixelMapY, pixelMap.pixelMapZ, p_x');
%     subplot(132); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, p_y');
%     subplot(133); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapY, p_z);
% end



%%
%% Create tubes

tube1 = PhantomTube(...
        [0, -0.2*mm, 70*mm], ...            % Center
        [200*um, 200*um, 30*mm], ...    % Size
        [0, 50, 90] ...                % Rotation
);
tube2 = PhantomTube(...
        [0, 0.2*mm, 70*mm], ...            % Center
        [200*um, 200*um, 30*mm], ...    % Size
        [0, -50, 90] ...                 % Rotation
);

% place bubbels
fps = 500;
dt = 1/fps;
dt = 0.05;
figure;
ha = axes(); xlabel('X'); ylabel('Y'); zlabel('Z'); hold(ha, 'on');
tube1.plot_skeleton(ha);
tube2.plot_skeleton(ha);

transducer.tx_aperture.plot_aperture(ha);
space.plot_skeleton(ha);

% pre-sim
tube1.update(1);
tube2.update(1);
view(3);
axis image

tube1.bub_scat.bub(6:end,:) = [];
tube1.bub_scat.pos(6:end,:) = [];
tube2.bub_scat.bub(6:end,:) = [];
tube2.bub_scat.pos(6:end,:) = [];

bub1_p = tube1.plot_bub_scat(ha);
bub2_p = tube2.plot_bub_scat(ha);

% pre-run
for fn = 1:UserSet.totalFrame
    tube1.update(dt);
    tube2.update(dt);
end


% Simulate
for f = 1:UserSet.totalFrame
    
%     waitbar(f/UserSet.totalFrame, ['Frame ',num2str(f),'/',num2str(UserSet.totalFrame)]);
    tube1.update(dt);
    tube2.update(dt);

    tube1.update_bub_scat(bub1_p);
    tube2.update_bub_scat(bub2_p);

    pause
end