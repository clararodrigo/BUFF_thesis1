function [img_out, pixelMap] = beamform_with_pxm(rf_data, coords, transducer, angles, focus)
    % create pixelmap
    pixelMap = mk_pxmap(coords);
    % beamform
    [img_out] = apply_cuDAS_better(rf_data, transducer, angles, focus, pixelMap);
end


function pxmap = mk_pxmap(coords)
    % buff crosstube
    pixelMap.dz=100e-6;% matrix
    pixelMap.dx=pixelMap.dz;% matrix
    pixelMap.dy=pixelMap.dx;
    pixelMap.upperLeft = coords(1,:);
    %     scatter(end,3)*1e3-5]/1000;%15 29 61
    pixelMap.bottomRight = coords(2,:);
    
    pixelMap.pixelMapX = pixelMap.upperLeft(1):pixelMap.dx:pixelMap.bottomRight(1);
    pixelMap.pixelMapY = pixelMap.upperLeft(2):pixelMap.dy:pixelMap.bottomRight(2);
    pixelMap.pixelMapZ = pixelMap.upperLeft(3):pixelMap.dz:pixelMap.bottomRight(3);
    
    pxmap = pixelMap;
    
    
end
