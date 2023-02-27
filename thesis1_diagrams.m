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

delays_x = repmat([-0.05:32:0.05],32,1);
plot3(delays_x,delays_x,reshape(delays,[32,32])*1000000)


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
title('Cross-tube setup')

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

%% Plot SRA
sra = load('saved100CompRndApod.mat');
sra = sra.CompRndApod(1,:,:);
sra1 = reshape(sra(1,1,:),[32,32]);
% sra2 = reshape(sra(1,2,:),[32,32]);
% sra3 = reshape(sra(1,3,:),[32,32]);
% sra4 = reshape(sra(1,4,:),[32,32]);

transducer.set_receive_apodization(1,1);

data = xdc_get(1, 'rect');
x = data([11, 20, 17, 14], :);
y = data([12, 21, 18, 15], :);
z = data([13, 22, 19, 16], :);
c = data([5, 5, 5, 5], :);

figure;
ha = axes(); xlabel('X'); ylabel('Y'); zlabel('Z'); hold(ha, 'on');
patch(ha, x, y, z, c);
view(ha, 3)
xlabel(ha, 'X');
ylabel(ha, 'Y');
zlabel(ha, 'Z');


