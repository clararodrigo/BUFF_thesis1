clear all;

addpath(genpath('RF'));
addpath(genpath('BF'));

%%

SDW_name = 'bf_0_SDW_*3cm.mat';
CROSS_name = 'bf_0_CROSS_*3cm.mat';

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

px = CROSS_data(4).pixelMap;

figure;
subplot(131); imagesc(px.pixelMapY, px.pixelMapZ, p_x');
subplot(132); imagesc(px.pixelMapX, px.pixelMapZ, p_y');
subplot(133); imagesc(px.pixelMapX, px.pixelMapY, p_z);

%% Normalize
for i = 1 : size(SDW_data,2)
    SDW_data(i).data = abs(SDW_data(i).data)./max(abs(SDW_data(i).data),[],'all');
    CROSS_data(i).data = abs(CROSS_data(i).data)./max(abs(CROSS_data(i).data),[],'all');
end

%% loop through and get FWHM
%  bubble 1 at z = 51 pixels ( 6cm depth )
%  bubble 2 at z = 151 pixels ( 7cm depth )
%  bubble 3 at z = 251 pixels ( 8cm depth )

for i = 1 : size(SDW_data,2)
    SDW_data(i).bub1 = SDW_data(i).data(:,:,51);        
    SDW_data(i).bub2 = SDW_data(i).data(:,:,151);
    SDW_data(i).bub3 = SDW_data(i).data(:,:,251);

    SDW_data(i).mid1 = sum(SDW_data(i).bub1, 2);            % calculate intensity profile of bub1
    SDW_data(i).mid2 = sum(SDW_data(i).bub2, 2);            % calculate intensity profile of bub1
    SDW_data(i).mid3 = sum(SDW_data(i).bub3, 2);            % calculate intensity profile of bub1

    % for CROSS
    CROSS_data(i).bub1 = CROSS_data(i).data(:,:,51);        
    CROSS_data(i).bub2 = CROSS_data(i).data(:,:,151);
    CROSS_data(i).bub3 = CROSS_data(i).data(:,:,251);

    CROSS_data(i).mid1 = sum(CROSS_data(i).bub1, 2);        % calculate intensity profile of bub1
    CROSS_data(i).mid2 = sum(CROSS_data(i).bub2, 2);        % calculate intensity profile of bub1
    CROSS_data(i).mid3 = sum(CROSS_data(i).bub3, 2);        % calculate intensity profile of bub1
    
end    

%%
plot_bubbles(SDW_data, 'SDW')
plot_bubbles(CROSS_data, 'CROSS')
%% find FWHM - SDW

SDW_FWHM = struct([]);
CROSS_FWHM = struct([]);
for i = 1 : size(SDW_data, 2)
    tmp_max = max(SDW_data(i).mid1,[],'all');
    tmp1 = find(fix(SDW_data(i).mid1) == round(tmp_max/2), 1, 'first');
    tmp2 = find(fix(SDW_data(i).mid1) == round(tmp_max/2), 1, 'last');
    SDW_FWHM(i).bub1 = [tmp1, tmp2];

    tmp_max = max(SDW_data(i).mid2,[],'all');
    tmp1 = find(fix(SDW_data(i).mid2) == round(tmp_max/2), 1, 'first');
    tmp2 = find(fix(SDW_data(i).mid2) == round(tmp_max/2), 1, 'last');
    SDW_FWHM(i).bub2 = [tmp1, tmp2];

    tmp_max = max(SDW_data(i).mid3,[],'all');
    tmp1 = find(fix(SDW_data(i).mid3) == round(tmp_max/2), 1, 'first');
    tmp2 = find(fix(SDW_data(i).mid3) == round(tmp_max/2), 1, 'last');
    SDW_FWHM(i).bub3 = [tmp1, tmp2];
end
% the i == 2 foesnt have a point exactly at half max
SDW_FWHM(2).bub1 = [185, 217];
SDW_FWHM(1).bub2 = [182, 220];
SDW_FWHM(3).bub1 = [186, 217];
SDW_FWHM(2).bub2 = [183, 219];

% convert to mm
for i = 1 : size(SDW_FWHM, 2)
    SDW_data(1,i).fwhm1 = px(SDW_FWHM(i).bub1(2)) - px(SDW_FWHM(i).bub1(1));
    SDW_data(1,i).fwhm2 = px(SDW_FWHM(i).bub2(2)) - px(SDW_FWHM(i).bub2(1));
    SDW_data(1,i).fwhm3 = px(SDW_FWHM(i).bub3(2)) - px(SDW_FWHM(i).bub3(1));
end

%% find FWHM - CROSS

CROSS_FWHM = struct([]);
for i = 1 : size(CROSS_data, 2)
    tmp_max = max(CROSS_data(i).mid1,[],'all');
    tmp1 = find(fix(CROSS_data(i).mid1) == round(tmp_max/2), 1, 'first');
    tmp2 = find(fix(CROSS_data(i).mid1) == round(tmp_max/2), 1, 'last');
    CROSS_FWHM(i).bub1 = [tmp1, tmp2];

    tmp_max = max(CROSS_data(i).mid2,[],'all');
    tmp1 = find(fix(CROSS_data(i).mid2) == round(tmp_max/2), 1, 'first');
    tmp2 = find(fix(CROSS_data(i).mid2) == round(tmp_max/2), 1, 'last');
    CROSS_FWHM(i).bub2 = [tmp1, tmp2];

    tmp_max = max(CROSS_data(i).mid3,[],'all');
    tmp1 = find(fix(CROSS_data(i).mid3) == round(tmp_max/2), 1, 'first');
    tmp2 = find(fix(CROSS_data(i).mid3) == round(tmp_max/2), 1, 'last');
    CROSS_FWHM(i).bub3 = [tmp1, tmp2];
end
% the i == 2 foesnt have a point exactly at half max
CROSS_FWHM(4).bub1 = [184, 218];
CROSS_FWHM(5).bub1 = [185, 217];
CROSS_FWHM(2).bub3 = [180, 222];
CROSS_FWHM(3).bub1 = [185, 217];
CROSS_FWHM(2).bub1 = [185, 217];

% convert to mm
for i = 1 : size(CROSS_FWHM, 2)
    CROSS_data(i).fwhm1 = px(CROSS_FWHM(i).bub1(2)) - px(CROSS_FWHM(i).bub1(1));
    CROSS_data(i).fwhm2 = px(CROSS_FWHM(i).bub2(2)) - px(CROSS_FWHM(i).bub2(1));
    CROSS_data(i).fwhm3 = px(CROSS_FWHM(i).bub3(2)) - px(CROSS_FWHM(i).bub3(1));
end
%% MLSLR
% From Bernal, 2020
% The MLSLR was calculated as the difference between the amplitude of 
% the main lobe (in dB) and the average background value

% getting it in dB - calculating ratio
thr = 0.2;
for i = 1 : size(SDW_data,2)
    tmp = SDW_data(i).bub1;
    tmp(tmp > thr) = 0;                                     % delete foreground pixels
    tmp_bckg = mean(tmp(~isnan(tmp)));                      % the mean of the background around bub at z=6.com
    SDW_data(i).db1 = SDW_data(i).bub1(200,200)/tmp_bckg;          


    tmp = SDW_data(i).bub2;
    tmp(tmp > thr) = 0;                                     % delete foreground pixels
    tmp_bckg = mean(tmp(~isnan(tmp)));                      % the mean of the background around bub at z=6.com
    SDW_data(i).db2 = SDW_data(i).bub2(200,200)/tmp_bckg;   

    tmp = SDW_data(i).bub3;
    tmp(tmp > thr) = 0;                                         % delete foreground pixels
    tmp_bckg = mean(tmp(~isnan(tmp)));                          % the mean of the background around bub at z=6.com
    SDW_data(i).db3 = SDW_data(i).bub3(200,200)/tmp_bckg;   
end

for i = 1 : size(CROSS_data,2)
    tmp = CROSS_data(i).bub1;
    tmp(tmp > thr) = 0;                                         % delete foreground pixels
    tmp_bckg = mean(tmp(~isnan(tmp)));                          % the mean of the background around bub at z=6.com
    CROSS_data(i).db1 = CROSS_data(i).bub1(200,200)/tmp_bckg;          


    tmp = CROSS_data(i).bub2;
    tmp(tmp > thr) = 0;                                         % delete foreground pixels
    tmp_bckg = mean(tmp(~isnan(tmp)));                          % the mean of the background around bub at z=6.com
    CROSS_data(i).db2 = CROSS_data(i).bub2(200,200)/tmp_bckg;   

    tmp = CROSS_data(i).bub3;
    tmp(tmp > thr) = 0;                                         % delete foreground pixels
    tmp_bckg = mean(tmp(~isnan(tmp)));                          % the mean of the background around bub at z=6.com
    CROSS_data(i).db3 = CROSS_data(i).bub3(200,200)/tmp_bckg;   
end


function plot_bubbles(data, name)
    figure;
    subplot(131)
    for i=1:5
        plot(data(i).pixelMap.pixelMapX, data(i).mid1); hold on
    end
    title('z = 6cm')
    legend('1 ap','2 ap','3 ap','4 ap','full ap');
    
    subplot(132)
    for i=1:5
        plot(data(i).pixelMap.pixelMapX, data(i).mid2); hold on
        ylim([0,50]);
    end
    title('z = 7cm')
    legend('1 ap','2 ap','3 ap','4 ap','full ap')
    
    subplot(133)
    for i=1:5
        plot(data(i).pixelMap.pixelMapX, data(i).mid3); hold on
        ylim([0,50]);
    end
    title('z = 8cm')
    legend('1 ap','2 ap','3 ap','4 ap','full ap')

    sgtitle(name);
end