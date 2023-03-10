%% Cleanup
clear all;
close all;
clc;
cd("/home/clararg/Documents/Scripts/Simulations/BUFF_thesis1")
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

UserSet.totalFrame = 200;
UserSet.ap = 2;
UserSet.old_bubbles = 1;


filename = ['CROSS_', num2str(UserSet.totalFrame), 'F', num2str(UserSet.ap), 'A_X'];

%% Setup
% Main volume
space = Box( ...
        [0, 0, 70*mm], ...              % Center
        [10*mm, 10*mm, 10*mm], ...      % Size
        [0, 0, 0] ...                   % Rotation
);


transducer = Matrix1024();
transducer.initial_load_apodizations('saved500CompRndApod.mat');
transducer.tx_aperture.excitation.n_cycles = 1;
transducer.tx_aperture.apply_waveforms();
transducer.tx_aperture.apply_delays();
transducer.tx_aperture.apply_apodization();
transducer.set_MI(0.05, space.center);

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
focus = -focus_r;

%% Create tubes

if(UserSet.old_bubbles)
    load('../RF/SDW/rf_SDW_200F2A_X_dt001_bubbles.mat')
else
    
    tube1 = PhantomTube(...
            [0, -0.2*mm, 70*mm], ...            % Center
            [200*um, 200*um, 15*mm], ...    % Size
            [0, 50, 90] ...                % Rotation
    );
    tube2 = PhantomTube(...
            [0, 0.2*mm, 70*mm], ...            % Center
            [200*um, 200*um, 15*mm], ...    % Size
            [0, -50, 90] ...                 % Rotation
    );
    
    % place bubbels
    % fps = 500;
    % dt = 1/fps;
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
    
    if(~UserSet.old_bubbles)
        tube1.update(dt);
        tube2.update(dt);
    
        tube1.update_bub_scat(bub1_p);
        tube2.update_bub_scat(bub2_p);
    
        tot_bub(:,f) = tube1.bub_scat + tube2.bub_scat;
    end
        
    fprintf('Frame %d/%d',f,UserSet.totalFrame);

    for a = 1 : size(focus_cart,1)
        fprintf('.');
        % steering
        delays = transducer.tx_aperture.calc_delays_diverging_pos(focus_cart(a,:));
        for ee = 1:transducer.tx_aperture.n_elements
            transducer.tx_aperture.elements(ee).delay = delays(ee);
        end
        transducer.tx_aperture.apply_delays();
        
        for ap = 1:UserSet.ap
            % configure transmission
            fprintf('.');
            transducer.set_transmit_apodization(f+a);
            transducer.set_receive_apodization(f+a, ap);
            
            bub_rf(f,a,ap) = transducer.tx_rx(tot_bub(:,f));
        end
    end
    if(mod(f,10)==0)
        tmp = bub_rf(f-9:f,:,:);
        save(['RF/bub_rf_CROSS_3A',num2str(f),'F'],'bub_rf','-v7.3');
    end
    fprintf('\n');
end

%% Load up data if deleted
if(~exist('bub_rf'))

    files = dir(['../RF/CROSS/',filename,'_*','.mat']);
    
    files(size(files,1)+1:size(files,1)+10) = files(2:11,:);
    files(2:11,:) = [];
    files(21,:) = files(3,:);
    files(3,:) = [];

    for i = 7 : 10                                       % I'll do it in thirds so my computer doesn't die
        load(['../RF/CROSS/', files(i).name]);
        if(i == 7)
            bub_rf = tmp;
        else
            bub_rf(size(bub_rf,1)+1 : size(bub_rf,1)+10,:,:) = tmp;
        end
    end
end
clear tmp;

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
    bub_rf(f,:,:) = zeros(10,2);
end
rf_data = permute(rf_data, [2,1,3,4,5]); % depth, elem, frames, angles, ap

% aperture compounding
rf_data = sum(rf_data,5);

%% to beamform in gpu
% save(['../RF/rf_', filename,'dt001'], 'rf_data', '-v7.3');
% save(['../RF/rf_', filename,'bubbles_dt001'], 'tot_bub','-v7.3');

%% beamforming
% save(['../RF/CROSS/rf_', filename,'_61_100'], 'rf_data','-v7.3');
load(['../RF/CROSS/rf_', filename,'_61_100'], 'rf_data');

%%
% ImgData = ones(301,301,301,30);
for f = 1:10:31 % 40:UserSet.totalFrame-39
    [tmp, pixelMap] = beamform_sim_ple(squeeze(rf_data(:,:,f:f+9,:)), transducer, focus_cart, focus);               % beamform
    ImgData(:,:,:,f:f+9) = permute(tmp, [3,2,1,4]);                                                                 % should just be 1 volume
end

%%
img = sum(ImgData,4);
img = permute(img, [3,2,1]);

img = abs(img);
p_x = squeeze(sum(img,1));
p_y = squeeze(sum(img,2));
p_z = squeeze(sum(img,3));

figure;
subplot(131); imagesc(pixelMap.pixelMapY, pixelMap.pixelMapZ, p_x');
subplot(132); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, p_y');
subplot(133); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapY, p_z);

%%
% saveIQ_simple(tmp, ['../BF/CROSS/bf_', filename,'_01_30'], pixelMap, UserSet)
data1 = getByteStreamFromArray(ImgData(:,:,:,1:15));
save(['../BF/CROSS/bf',filename,'61_75'],'data1','-v7.3')

data2 = getByteStreamFromArray(ImgData(:,:,:,16:30));
save(['../BF/CROSS/bf',filename,'76_90'],'data2','-v7.3')

data3 = getByteStreamFromArray(ImgData(:,:,:,31:40));
save(['../BF/CROSS/bf',filename,'90_100'],'data3','-v7.3')
