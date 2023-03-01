function [img_out, pixelMap] = beamform_sim_ple(rf_data, transducer, angles, focus)
    % create pixelmap
    pixelMap = mk_pxmap();
    % beamform
    [img_out] = apply_cuDAS_better(rf_data, transducer, angles, focus, pixelMap);
end


function pxmap = mk_pxmap()
    % buff crosstube
    pixelMap.dz=100e-6;% matrix
    pixelMap.dx=pixelMap.dz;% matrix
    pixelMap.dy=pixelMap.dx;
    pixelMap.upperLeft = [
        -15e-3, ...                                   % 15 for X, 20 for contrast phantom                  
        -15e-3, ...
         55e-3 ...                                    % 50-70 for contrast, 55-85 for X
    ];
    %     scatter(end,3)*1e3-5]/1000;%15 29 61
    pixelMap.bottomRight = [...
        15e-3,...
        15e-3,...
        85e-3];
    
    pixelMap.pixelMapX = pixelMap.upperLeft(1):pixelMap.dx:pixelMap.bottomRight(1);
    pixelMap.pixelMapY = pixelMap.upperLeft(2):pixelMap.dy:pixelMap.bottomRight(2);
    pixelMap.pixelMapZ = pixelMap.upperLeft(3):pixelMap.dz:pixelMap.bottomRight(3);
    
    pxmap = pixelMap;
    
    
end
