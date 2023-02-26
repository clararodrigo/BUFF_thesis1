%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
addpath(genpath('D:\Share to student\01 Save Fast'));
%BUFF
addpath(genpath('F:\Clara\Simulations\Field\buff-0.1.0-alpha\src\'));
%FieldII
addpath(genpath("F:\Clara\Simulations\Field\m_files"));
addpath(genpath("F:\Clara\Simulations\Field\"));
%save directories
addpath(genpath('F:\Clara\Simulations\Field\matrix_array_sim\matrix_array_sim\RF_SDW'));
addpath(genpath('F:\Clara\Simulations\Field\matrix_array_sim\matrix_array_sim\BF_SDW'));

field_init(0);
set_sampling(GlobalConfig().fs);

UserSet.totalFrame = 1;
UserSet.ap = 2;
filename = ['CROSS_1F', num2str(UserSet.ap), 'A_3bub'];

%% Setup
% Main volume
space = Box( ...
        [0, 0, 50*mm], ...              % Center
        [10*mm, 10*mm, 10*mm], ...      % Size
        [0, 0, 0] ...                   % Rotation
);


transducer = Matrix1024();
transducer.initial_load_apodizations('F:\Clara\Simulations\Field\saved100CompRndApod.mat');
transducer.tx_aperture.excitation.n_cycles = 1;
transducer.tx_aperture.apply_waveforms();
transducer.tx_aperture.apply_delays();
transducer.tx_aperture.apply_apodization();
transducer.set_MI(0.05, space.center);

%% Creeate 3 bubbles
bub = SonovueBubble(3.6e-6);
bub_scat = BubbleScatterers([space.center; space.center+[0,0,10*mm];...
    space.center+[0,0,20*mm]], [bub,bub,bub]);
% bub_scat = BubbleScatterers([space.center], [bub]);
%% Angle information
lambda = transducer.c/(transducer.f0);


x_focus = (-50:25:50)*(32/150);
y_focus = (-50:25:50)*(32/150);
z_focus = 32:32;

% % Grid
%[XF, YF, ZF] = ndgrid(x_focus, y_focus, z_focus);

% Cross
XF = [x_focus  , x_focus *0 ] * lambda;
YF = [y_focus *0  , y_focus ] * lambda;
ZF = ones(size(XF)) *32 * lambda;

% Cartesian
focus_cart = [XF(:), YF(:), ZF(:)];

% Spherical
focus_r = vecnorm(focus_cart,2,2);
focus_a2 = asin(YF(:)./focus_r);
focus_a1 = asin(XF(:)./focus_r./cos(focus_a2));
angles = [focus_a1(:), focus_a2(:)];
focus = -focus_r;

%%
% ha = axes();
% hold all;
% transducer.tx_aperture.plot_aperture(ha);
% plot3(ha, XF(:), YF(:), ZF(:), 'ro');


for f = 1:UserSet.totalFrame
    fprintf('Frame %d/%d :',f,UserSet.totalFrame);
    for a = 1:size(focus_cart,1)
        fprintf('|');   
        
        % steering
        delays = transducer.tx_aperture.calc_delays_diverging_pos(focus_cart(a,:));
%         delays = transducer.tx_aperture.calc_delays(angles(a,:));
        for ee = 1:transducer.tx_aperture.n_elements
            transducer.tx_aperture.elements(ee).delay = delays(ee);
        end
        transducer.tx_aperture.apply_delays();
            
        
        for ap = 1:UserSet.ap            
            fprintf('.');
            % apply apodizations
            transducer.set_transmit_apodization(f+a);
            transducer.set_receive_apodization(f+a, ap);

            bub_rf(f,a,ap) = transducer.tx_rx(bub_scat);
        end
    end
    fprintf('\n');
end

%% Diagram
figure();
ha = axes();
hold(ha, 'on');
transducer.tx_aperture.plot_aperture(ha);
space.plot_skeleton(ha)
% Scatterers
% plot3(bub_scat.pos(:,1),bub_scat.pos(:,2),bub_scat.pos(:,3),'ob', 'MarkerSize', 20);

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

% angle compounding
rf_data = sum(rf_data,5);
% rf_data = sum(rf_data,3);

%% to beamform in gpu
% save(['PSF/rf_', filename], 'rf_data');
ImgData = zeros(201,201,301,size(angles,1));
for a = 1 : size(angles,1)
    [tmp, pixelMap] = beamform_sim_ple(rf_data(:,:,1,1), transducer, angles, focus);
    ImgData(:,:,:,a) = permute(tmp, [3,2,1]);
end
%%
% save data
saveIQ_simple(ImgData, ['PSF/bf_', filename], pixelMap, UserSet)
%%
img = sum(ImgData,4);
% img = permute(img, [3,2,1]);

img = abs(img);
p_x = squeeze(sum(img,1));
p_y = squeeze(sum(img,2));
p_z = squeeze(sum(img,3));

figure;
subplot(131); imagesc(pixelMap.pixelMapY, pixelMap.pixelMapZ, p_x');
subplot(132); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, p_y');
subplot(133); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapY, p_z);

%%
t = permute(img_out, [3,2,1,4,5]);
img = abs(t(:,:,:,1,1));
p_x = squeeze(sum(img,1));
p_y = squeeze(sum(img,2));
p_z = squeeze(sum(img,3));

figure;
subplot(131); imagesc(pixelMap.pixelMapY, pixelMap.pixelMapZ, p_x');
subplot(132); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapZ, p_y');
subplot(133); imagesc(pixelMap.pixelMapX, pixelMap.pixelMapY, p_z);

%%
img = abs(img)./max(abs(img),[],'all');
mid = sum(img(:,:,1),2);
