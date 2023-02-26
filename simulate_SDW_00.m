%% Cleanup
clear all;
close all;
clc;

%% Import Libraries
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
UserSet.ap = 4;

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

%% Creeate one bubble
bub = SonovueBubble(3.6e-6);
bub_scat = BubbleScatterers(space.center, [bub]);

%% Simulate

delays = transducer.tx_aperture.calc_delays_diverging_pos([0,0,-0.01]);
for ee = 1:transducer.tx_aperture.n_elements
    transducer.tx_aperture.elements(ee).delay = delays(ee);
end
transducer.tx_aperture.apply_delays();

for ap = 1:UserSet.ap
    % configure transmission
    transducer.set_transmit_apodization(1);
    transducer.set_receive_apodization(1, ap);
    bub_rf(ap) = transducer.tx_rx(bub_scat);
end


%% Diagram
figure();
ha = axes();
hold(ha, 'on');
transducer.tx_aperture.plot_aperture(ha);
space.plot_skeleton(ha)
% Scatterers
plot3(bub_scat.pos(:,1),bub_scat.pos(:,2),bub_scat.pos(:,3),'ob', 'MarkerSize', 20);

%% Pad signals 

max_t = max(arrayfun(@(x) max(bub_rf(x).time_vector), 1:3));

for si = 1:numel(bub_rf)
    bub_rf(si) = bub_rf(si) + RFSignals(zeros(1024,2), 0) + RFSignals(zeros(1024,2), max_t);
end

%% Cat RF data for all apodizations
rf_data = cat(3,bub_rf.signals);
rf_data = permute(rf_data, [2,1,3]);

% compounding
rf_data = sum(rf_data,3);

%% to beamform in gpu
% save('../RF_SDW/rf_1F4A_3bub','rf_data');
myMatrix_Pre_BF_clara;
%%
figure;
imagesc(image_coords(:,1),image_coords(:,2),squeeze(abs(scat_bf(:,:))))

%%
img = abs(img)./max(abs(img),[],'all');
mid = sum(img(:,:,1),2);
