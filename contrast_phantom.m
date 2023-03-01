%  By Clara: creating a phantom to assess CNR
%  will place high-scatterers at x = -1cm, y = 0cm
%  and low-scatterers at x = 1cm, y = 0cm

function [positions, amp] = contrast_phantom (N)
    x_size = 40/1000;   %  Width of phantom [mm]
    y_size = 10/1000;   %  Transverse width of phantom [mm]
    z_size = 20/1000;   %  Height of phantom [mm]
    z_start = 50/1000;  %  Start of phantom surface [mm];
    
    %  Create the general scatterers
    
    x = (rand (N,1)-0.5)*x_size;
    y = (rand (N,1)-0.5)*y_size;
    z = rand (N,1)*z_size + z_start;
    
    %%  Generate the amplitudes with a Gaussian distribution
    
%     amp = randn(N,1)*2;
    amp = ones(N,1);
    
    %  Make the cyst and set the amplitudes to zero inside
    
    r = 5*mm;                                        %  Radius of cyst [mm]

    xc = -15*mm;                                     
    zc = 10/1000 + z_start;    
    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside) + 0.2*amp .* inside;

    xc = -5*mm;                                     
    zc = 10/1000 + z_start;    
    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside) + 0.3*amp .* inside; 

    xc = 5*mm;                                     
    zc = 10/1000 + z_start;    
    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside) + 0.5*amp .* inside; 

    xc = 15*mm;                                     
    zc = 10/1000 + z_start;    
    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside) + 0.7*amp .* inside; 


    %  Return the variables
    positions=[x y z];
end