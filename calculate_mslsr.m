clear all;

addpath(genpath('RF'));
addpath(genpath('BF'));

%%
nor = 0;                                                    % use this to choose whether well normalize stuff 

SDW_name = 'bf_2_SDW_*A.mat';
CROSS_name = 'bf_2_CROSS_*A.mat';

SDW_files = dir(['../BF/SDW/',SDW_name]);
CROSS_files = dir(['../BF/CROSS/',CROSS_name]);

%%
SDW_data = struct([]);
CROSS_data = struct([]);
for i = 1 : size(SDW_files,1)
    tmp = load(['../BF/SDW/', SDW_files(i).name]);
    SDW_data(i).data = getArrayFromByteStream(tmp.data);
    SDW_data(i).pixelMap = tmp.pixelMap;
    SDW_data(i).UserSet = tmp.UserSet;
    
    tmp = load(['../BF/CROSS/', CROSS_files(i).name]);
    CROSS_data(i).data = getArrayFromByteStream(tmp.data);
    CROSS_data(i).pixelMap = tmp.pixelMap;
    CROSS_data(i).UserSet = tmp.UserSet;
end

%% plotting one to see
img_tmp = abs(CROSS_data(1).data);
p_x = squeeze(sum(img_tmp,1));
p_y = squeeze(sum(img_tmp,2));
p_z = squeeze(sum(abs(img_tmp),3));

px = CROSS_data(1).pixelMap;

figure;
subplot(131); imagesc(px.pixelMapY, px.pixelMapZ, p_x');
subplot(132); imagesc(px.pixelMapX, px.pixelMapZ, p_y');
subplot(133); imagesc(px.pixelMapX, px.pixelMapY, p_z);

%% Normalize
if nor
    for i = 1 : size(SDW_data,2)
        SDW_data(i).data = abs(SDW_data(i).data)./max(abs(SDW_data(i).data),[],'all');
        CROSS_data(i).data = abs(CROSS_data(i).data)./max(abs(CROSS_data(i).data),[],'all');
    end
else
    for i = 1 : size(SDW_data,2)
        SDW_data(i).data = abs(SDW_data(i).data);
        CROSS_data(i).data = abs(CROSS_data(i).data);
    end
end
%% loop through and get FWHM

loc = 11 : 10 : 501;                                                        % these are the z-locations of the bubbles

for i = 1 : size(SDW_data,2)
    for bub = 1 : 50
        bub_name = ['bub',num2str(bub)];
        mid_name = ['mid',num2str(bub)];

        SDW_data(i).(bub_name) = SDW_data(i).data(:,:,loc(bub));   
        SDW_data(i).(mid_name) = sum(SDW_data(i).(bub_name), 2);            % calculate intensity profile 

        % for CROSS
        CROSS_data(i).(bub_name) = CROSS_data(i).data(:,:,loc(bub));        
        CROSS_data(i).(mid_name) = sum(CROSS_data(i).(bub_name), 2);        % calculate intensity profile of bub1
    end
end    

%%
plot_bubbles(SDW_data, 'SDW')
plot_bubbles(CROSS_data, 'CROSS')

%% MLSLR
% From Bernal, 2020
% The MLSLR was calculated as the difference between the amplitude of 
% the main lobe (in dB) and the average background value

% getting it in dB - calculating ratio
thr = 0.2;
for i = 1 : size(SDW_data,2)
    for bub = 1 : 50


        bub_name = ['bub',num2str(bub)];
        tmp = SDW_data(i).(bub_name);
        tmp(tmp > thr) = 0;                                                 % delete foreground pixels
        tmp_bckg = mean(tmp(~isnan(tmp)));                                  % the mean of the background around bub at z=6.com
        SDW_data(i).(['db',num2str(bub)]) = SDW_data(i).(bub_name)(200,200)/tmp_bckg; 

        bub_name = ['bub',num2str(bub)];
        tmp = CROSS_data(i).(bub_name);
        tmp(tmp > thr) = 0;                                                 % delete foreground pixels
        tmp_bckg = mean(tmp(~isnan(tmp)));                                  % the mean of the background around bub at z=6.com
        CROSS_data(i).(['db',num2str(bub)]) = CROSS_data(i).(bub_name)(200,200)/tmp_bckg; 
    end
end

%% Visualize MLSLR

val = ones(2,5,50);         % [ SDW/CROSS, #ap, bubble ]


for bub = 1 : 50
    name = ['db',num2str(bub)];
    val(1, :, bub) = [SDW_data.(name)];
    val(2, :, bub) = [CROSS_data.(name)];
end

figure;
subplot(121); 
for ap = 1 : 5
    semilogy([1:50],squeeze(val(1,ap,:))); hold on
%     plot(squeeze(val(1,ap,:))); hold on 
end
legend('1 ap','2 ap','3 ap','4 ap','full ap'); 
ylabel('MLSLR')
title('SDW');
% set(gca,'XTick',[0, 15, 30, 45, 60])
% set(gca,'XTickLabel',[3:1:8])

subplot(122); 
for ap = 1 : 5
    semilogy([1:50],squeeze(val(2,ap,:))); hold on
end
legend('1 ap','2 ap','3 ap','4 ap','full ap');
ylabel('MLSLR')
title('CROSS')

%%
function plot_bubbles(data, name)

    figure('Renderer', 'painters', 'Position', [10 10 900 800])
    subplot(131)
    for i=1:5
        plot(data(i).pixelMap.pixelMapX*100, data(i).mid1); hold on
    end
    title('z = 3cm')
    legend('1 ap','2 ap','3 ap','4 ap','full ap');
    
    subplot(132)
    for i=1:5
        plot(data(i).pixelMap.pixelMapX*100, data(i).mid20); hold on
%         ylim([0,50]);
    end
    title('z = 5cm')
    legend('1 ap','2 ap','3 ap','4 ap','full ap')
    
    subplot(133)
    for i=1:5
        plot(data(i).pixelMap.pixelMapX*100, data(i).mid50); hold on
%         ylim([0,50]);
    end
    title('z = 8cm')
    legend('1 ap','2 ap','3 ap','4 ap','full ap')

    sgtitle(name);
end