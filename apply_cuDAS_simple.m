function [img_out] = apply_cuDAS_simple(rf_data, transducer, angles, focus, pixelMap)

    addpath(genpath("D:\Share to student\Verasonics_GPU_Beamformer_CheeHau\src\gpuDAS"));

    
    nframes = size(rf_data, 3);
    ntrans = size(rf_data, 4);
    
    % Delay due to transducer impulse response
    delay = length(transducer.tx_aperture.impulse_response.signal)/transducer.fs;

    % Define pre-filtering FIR coefficients
    prefilt_fir = 1;

    % pre-allocate memory for output
    out_sz = [...
                length(pixelMap.pixelMapZ), ...
                length(pixelMap.pixelMapY), ...
                length(pixelMap.pixelMapX) ...
            ];
    img_out = zeros([out_sz, nframes, ntrans]);
    
    % Convert data type to simple
    if ~isa(rf_data, 'single')
        rf_data = single(rf_data);
    end
    
    focus = single(focus);
    angles = single(angles);
    
    for f = 1 : nframes
        disp(['Processing Plane Wave Frame #',num2str(f)]);
        for a = 1 : ntrans
            
            % Initialize GPU beamforming with parameters
            cuDAS_single(1, ...
                squeeze(rf_data(:,:,f,a)), ...                                    % data: [depth, elements, ...]
                angles(a, 1), ...                               % steering angles asimuth: [nfiringframes, 1]
                angles(a, 2), ...                               % steering angles elevation: [nfiringframes, 1]
                transducer.tx_aperture.elem_centers(:, 1), ...  % element centers X
                transducer.tx_aperture.elem_centers(:, 2), ...  % element centers Y
                transducer.tx_aperture.elem_centers(:, 3), ...  % element centers Z
                prefilt_fir, ...                                % FIR filter coefficients pre beamforming
                pixelMap.pixelMapX, ...                         % pixelmap X
                pixelMap.pixelMapY, ...                         % pixelmap Y
                pixelMap.pixelMapZ, ...                         % pixelmap Z
                delay, ...                                      % Delay due to transducer impulse response
                transducer.fs, ...                              % Sampling frequency
                transducer.tx_aperture.f0, ...                  % Transmitted central frequency
                transducer.c, ...                               % speed of sound
                focus(a), ...                                      % bad comment [1,3] position for diverging wave delay calculation
                0, ...                                          % Sensitivity cutoff = 0 -> no cutoff
                0, ...                                          % 0: High Resolution Image
                0, ...                                          % GPUID 0 (first and only card)          
                2, ...                                          % 2: 3D
                0 ...                                           % Cartesian coordinates
            );

    
            disp(['Angle ',num2str(a)]);
            % Run Beamforming
            img_out(:, :, :, f, a) = cuDAS_single(0, squeeze(rf_data(:, :, f, a)));
            
            % de-initialize GPU
            cuDAS_single(-1);
        end
    end

    


end