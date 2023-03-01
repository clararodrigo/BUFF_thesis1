%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('../../Share to student/Verasonics_GPU_Beamformer_CheeHau/src/gpuDAS'));
addpath(genpath('/home/clararg/Documents/Scripts/Share to student/01 Save Fast'));
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
UserSet.full = 0;                                                               % if 1: not using sparse apertures
medium = 0;                                                                     % 0 : three single bubbles; 
                                                                                % 1: field ii contrast phantom

filename = [num2str(medium), '_CROSS_', num2str(UserSet.totalFrame), ...        % used for saving
            'F', num2str(UserSet.ap), 'A_3cm'];
if(UserSet.full)
    UserSet.ap = 1;
    filename = [num2str(medium), '_CROSS_', num2str(UserSet.totalFrame), ...    % used for saving
            'F5A_3cm'];
end

%% Setup
% Main volume
space = Box( ...
        [0, 0, 30*mm], ...              % Center
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
                        space.center+[0, 0*mm, 10*mm];...
                        space.center+[0, 0*mm, 20*mm]...
                        ], [bub,bub,bub]);
    case 1
        [pos, amp] = contrast_phantom(1000);
        scat = LinearScatterers(pos,amp);
end

%% Angle information
lambda = transducer.c/(transducer.f0);


x_focus = (-50:25:50)*(32/150);
y_focus = (-50:25:50)*(32/150);
z_focus = 32:32;

% Cross
XF = [x_focus  , x_focus *0 ] * lambda;
YF = [y_focus *0  , y_focus ] * lambda;
ZF = ones(size(XF)) *32 * lambda;

% Cartesian
focus_cart = [XF(:), YF(:), ZF(:)];

% Spherical
focus_r = vecnorm(focus_cart,2,2);
% focus_a2 = asin(YF(:)./focus_r);
% focus_a1 = asin(XF(:)./focus_r./cos(focus_a2));
% angles = [focus_a1(:), focus_a2(:)];
% angles = focus_cart;
focus = -focus_r;

%% Simulate

for f = 1:UserSet.totalFrame
    fprintf('Frame %d/%d :',f,UserSet.totalFrame);
    for a = 1:size(focus_cart,1)
        fprintf('|');   
        
        % steering
        delays = transducer.tx_aperture.calc_delays_diverging_pos(focus_cart(a,:));

        for ee = 1:transducer.tx_aperture.n_elements
            transducer.tx_aperture.elements(ee).delay = delays(ee);
        end
        transducer.tx_aperture.apply_delays();
            
        
        if(~UserSet.full), transducer.set_transmit_apodization(f+a); end    % if using SRAC, apply transmit apod
        for ap = 1:UserSet.ap            
            fprintf('.');
            
            if(~UserSet.full)
                transducer.set_receive_apodization(f+a, ap);                % if using SRAC, apply receive apod
            end

            bub_rf(f,a,ap) = transducer.tx_rx(scat);
        end
    end
    fprintf('\n');
end


%% Pad signals 
max_t = max(arrayfun(@(x) max(bub_rf(x).time_vector), 1:numel(bub_rf)));

for ii = 1:numel(bub_rf)
	bub_rf(ii) = RFSignals(zeros(1024,2), 0) + bub_rf(ii) + RFSignals(zeros(1024,2), max_t);
end

%% Cat RF data for all apodizations
for f = 1:UserSet.totalFrame
    for a = 1:10
        for ap = 1:UserSet.ap
            rf_data(:,:,f,a,ap) = bub_rf(f,a,ap).signals;
        end
    end
end
rf_data = permute(rf_data, [2,1,3,4,5]); % depth, elem, frames, angles, ap

% compounding
rf_data = sum(rf_data,5);

%% to beamform in gpu
save(['../RF/CROSS/rf_', filename], 'rf_data');

%%
load('../RF/CROSS/rf_0_CROSS_1F5A.mat')
filename = '0_CROSS_1F5A';
UserSet.ap=1;
UserSet.full=1;
% 
rf_data = squeeze(rf_data);                                                     % to match needed shape
[tmp, pixelMap] = beamform_sim_ple(rf_data, transducer, focus_cart, focus);     % beamform
ImgData = permute(tmp, [3,2,1,4]);                                              % should just be 1 volume

%%
img = sum(ImgData,4);

img = abs(img);
p_x = squeeze(sum(img,1));
p_y = squeeze(sum(img,2));
p_z = squeeze(sum(img,3));

figure;
subplot(131); imagesc(pixelMap.pixelMapY, pixelMap.pixelMapZ, p_x');
subplot(132); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, p_y');
subplot(133); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapY, p_z);

%% Save data
saveIQ_simple(ImgData, ['../BF/CROSS/bf_', filename], pixelMap, UserSet)