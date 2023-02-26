%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('../../Share to student/Verasonics_GPU_Beamformer_CheeHau/src/gpuDAS'));
%BUFF
addpath(genpath('../buff/src'));
%FieldII
addpath(genpath("../field_ii"));
%save directories
addpath(genpath('RF'));
addpath(genpath('BF'));

field_init(0);
set_sampling(GlobalConfig().fs);

UserSet.totalFrame = 1;
UserSet.ap = 4;
medium = 1;                                                             % 0 : three single bubbles; 
                                                                        % 1: field ii contrast phantom

filename = [num2str(medium), '_SDW_', num2str(UserSet.totalFrame), ...  % used for saving
            'F', num2str(UserSet.ap), 'A_'];

%% Setup
% Main volume
space = Box( ...
        [0, 0, 52.5*mm], ...              % Center
        [10*mm, 10*mm, 10*mm], ...      % Size
        [0, 0, 0] ...                   % Rotation
);


transducer = Matrix1024();
transducer.initial_load_apodizations('saved100CompRndApod.mat');
transducer.tx_aperture.excitation.n_cycles = 1;
transducer.tx_aperture.apply_waveforms();
transducer.tx_aperture.apply_delays();
transducer.tx_aperture.apply_apodization();
transducer.set_MI(0.05, space.center);

%% Create Medium
switch medium 
    case 0
        bub = SonovueBubble(3.6e-6);
        scat = BubbleScatterers([...
                        space.center; ...
                        space.center+[0, 10*mm, 10*mm];...
                        space.center+[0, 0,     20*mm]...
                        ], [bub,bub,bub]);
    case 1
        [pos, amp] = cyst_phantom(1000);
        scat = LinearScatterers(pos,amp);
end


%% Set focus and steering
lambda = transducer.c/(transducer.f0);

delays = transducer.tx_aperture.calc_delays_diverging_pos([0,0,-32*lambda]);
for ee = 1:transducer.tx_aperture.n_elements
    transducer.tx_aperture.elements(ee).delay = delays(ee);
end
transducer.tx_aperture.apply_delays();


%% Simulate
for f = 1:UserSet.totalFrame
    fprintf('Frame %d/%d',f,UserSet.totalFrame);
    
    transducer.set_transmit_apodization(f);             % apply transmit sparse apertures
    for ap = 1:UserSet.ap
        fprintf('.');

        transducer.set_receive_apodization(f, ap);      % apply receive sparse aperture
        bub_rf(f,ap) = transducer.tx_rx(scat);          % get rf signals
    end
    fprintf('\n');
end

%% Pad signals 

max_t = max(arrayfun(@(x) max(bub_rf(x).time_vector), 1:UserSet.ap));

for si = 1:numel(bub_rf)
    bub_rf(si) = bub_rf(si) + RFSignals(zeros(1024,2), 0) + RFSignals(zeros(1024,2), max_t);
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
% save(['../RF/SDW/rf_', filename], 'rf_data');
%%
lambda = transducer.c/(transducer.f0);
angles = [0,0];
focus = -32*lambda;
[ImgData, pixelMap] = beamform_sim_ple(rf_data, transducer, angles, focus);

% beamform yields ImgData in ZYX format, transform to XYZ
ImgData = permute(ImgData, [3,2,1]);

%%
% saveIQData(ImgData,'BF/','bf_SDW_50F2A_X_2',pixelMap,UserSet,Trans,ImgParam);

%%
ImgData = abs(ImgData);
p_x = squeeze(sum(ImgData,1));
p_y = squeeze(sum(ImgData,2));
p_z = squeeze(sum(abs(ImgData),3));

figure;
subplot(131); imagesc(pixelMap.pixelMapY, pixelMap.pixelMapZ, p_x');
subplot(132); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, p_y');
subplot(133); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapY, p_z);
%%
figure;
for i = 1 : size(ImgData,2)
    imagesc(abs(squeeze(ImgData(:,i,:)))');
    pause;
end