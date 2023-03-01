%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/Verasonics_GPU_Beamformer_CheeHau/src/gpuDAS'));
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/01 Save Fast'));
%BUFF
addpath(genpath('../buff/src'));
%FieldII
addpath(genpath("../field_ii"));
%save directories
addpath(genpath('../RF'));
addpath(genpath('../BF'));

field_init(0);
set_sampling(GlobalConfig().fs);

UserSet.totalFrame = 200;
UserSet.ap = 3;
UserSet.old_bubbles = 1;


filename = ['SDW_', num2str(UserSet.totalFrame), 'F', num2str(UserSet.ap), 'A_X'];

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
%% Create tubes
if(UserSet.old_bubbles)
    fprintf('Using previously defined tracks \n')
    load('../RF/SDW/rf_SDW_200F2A_X_dt001_bubbles.mat')
else
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
end

% Simulate
for f = 1:UserSet.totalFrame
    
%     waitbar(f/UserSet.totalFrame, ['Frame ',num2str(f),'/',num2str(UserSet.totalFrame)]);
    if(~UserSet.old_bubbles)
        tube1.update(dt);
        tube2.update(dt);
    
        tube1.update_bub_scat(bub1_p);
        tube2.update_bub_scat(bub2_p);
    
        tot_bub(:,f) = tube1.bub_scat + tube2.bub_scat;
    end
        
    fprintf('Frame %d/%d',f,UserSet.totalFrame);
    for ap = 1:UserSet.ap
        % configure transmission
        fprintf('.');
        transducer.set_transmit_apodization(f);
        transducer.set_receive_apodization(f, ap);
        bub_rf(f,ap) = transducer.tx_rx(tot_bub(:,f));
    end
    fprintf('\n');

    if(mod(f,10)==0)
        tmp = bub_rf(f-9:f,:,:);
        save(['RF/bub_rf_SDW_1A',num2str(f),'F'],'bub_rf','-v7.3');
        bub_rf(f,:) = zeros(1,2);
    end
    fprintf('\n');
end

%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/Verasonics_GPU_Beamformer_CheeHau/src/gpuDAS'));
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/01 Save Fast'));
%BUFF
addpath(genpath('../buff/src'));
%FieldII
addpath(genpath("../field_ii"));
%save directories
addpath(genpath('../RF'));
addpath(genpath('../BF'));

field_init(0);
set_sampling(GlobalConfig().fs);

UserSet.totalFrame = 200;
UserSet.ap = 1;
UserSet.old_bubbles = 0;


filename = ['SDW_', num2str(UserSet.totalFrame), 'F', num2str(UserSet.ap), 'A_X'];

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
%% Create tubes
if(UserSet.old_bubbles)
    fprintf('Using previously defined tracks \n')
    load('../RF/SDW/rf_SDW_200F2A_X_dt001_bubbles.mat')
else
    tube1 = PhantomTube(...
            [0, -0.2*mm, 70*mm], ...            % Center
            [200*um, 200*um, 30*mm], ...    % Size
            [0, 20, 90] ...                % Rotation
    );
    tube2 = PhantomTube(...
            [0, 0.2*mm, 70*mm], ...            % Center
            [200*um, 200*um, 30*mm], ...    % Size
            [0, -20, 90] ...                 % Rotation
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
end

% Simulate
for f = 1:UserSet.totalFrame
    
%     waitbar(f/UserSet.totalFrame, ['Frame ',num2str(f),'/',num2str(UserSet.totalFrame)]);
    if(~UserSet.old_bubbles)
        tube1.update(dt);
        tube2.update(dt);
    
        tube1.update_bub_scat(bub1_p);
        tube2.update_bub_scat(bub2_p);
    
        tot_bub(:,f) = tube1.bub_scat + tube2.bub_scat;
    end
        
    fprintf('Frame %d/%d',f,UserSet.totalFrame);
    for ap = 1:UserSet.ap
        % configure transmission
        fprintf('.');
        transducer.set_transmit_apodization(f);
        transducer.set_receive_apodization(f, ap);
        bub_rf(f,ap) = transducer.tx_rx(tot_bub(:,f));
    end
    fprintf('\n');

    if(mod(f,10)==0)
        tmp = bub_rf(f-9:f,:,:);
        save(['RF/bub_rf_SDW_1A',num2str(f),'F'],'bub_rf','-v7.3');
        bub_rf(f,:) = zeros(1,2);
    end
    fprintf('\n');
end