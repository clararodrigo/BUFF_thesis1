%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('../../Share to student/Verasonics_GPU_Beamformer_CheeHau/src/gpuDAS'));
addpath(genpath('../../Share to student/01 Save Fast'));
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
UserSet.ap = 2;
UserSet.full = 0;                                                           % if 1: not using sparse apertures
medium = 1;                                                                 % 0 : three single bubbles; 
                                                                            % 1: field ii contrast phantom

if(UserSet.full)
    UserSet.ap = 1;
    filename = [num2str(medium), '_SDW_', num2str(UserSet.totalFrame), ...  % used for saving
                'F', '5A'];
else
    filename = [num2str(medium), '_SDW_', num2str(UserSet.totalFrame), ...  % used for saving
            'F', num2str(UserSet.ap), 'A'];
end

%% Setup
% Main volume
space = Box( ...
        [0, 0, 60*mm], ...              % Center
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
                        space.center; ...                       % Location of bubble 1
                        space.center+[0, 0*mm, 10*mm];...       % Location of bubble 2
                        space.center+[0, 0*mm, 20*mm]...        % Location of bubble 3
                        ], [bub, bub, bub]);
    case 1
        [pos, amp] = contrast_phantom(10000);
        scat = LinearScatterers(pos,amp);

    case 2
        bub = SonovueBubble(3.6e-6);
        
        positions = [];
        for i = 1:50
            positions = [positions; space.center+[0, 0, i*1*mm]];

        end

        scat = BubbleScatterers(positions, repmat([bub],1,50));
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
    
    if(~UserSet.full), transducer.set_transmit_apodization(f); end    % if using SRAC, apply transmit apod
    for ap = 1:UserSet.ap
        fprintf('.');

        if(~UserSet.full)
            transducer.set_transmit_apodization(f);
                transducer.set_receive_apodization(f, ap);                % if using SRAC, apply receive apod
        end
        
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

%% Beamform images
% load('RF/SDW/rf_0_SDW_1F5A_3cm.mat')
% filename = '0_SDW_1F5A_3cm';
% UserSet.ap = 1;
% UserSet.full = 1;

lambda = transducer.c/(transducer.f0);
angles = [0,0];
focus = -32*lambda;
coords = [[-25*mm, -5*mm, 50*mm];[25* mm, 5*mm, 70*mm]];
[ImgData, pixelMap] = beamform_with_px(rf_data, coords, transducer, angles, focus);

% beamform yields ImgData in ZYX format, transform to XYZ
ImgData = permute(ImgData, [3,2,1]);


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
% saveIQ_simple(ImgData,['../BF/SDW/bf_', filename],pixelMap,UserSet); 
figure;
for i = 1:10:200 
    imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, squeeze(abs(ImgData(i,:,:)))); pause;
end
