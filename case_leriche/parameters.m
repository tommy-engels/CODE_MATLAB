function parameters
    global params
    params.CASE = 'leriche';
    params.Lx   = 128*2.0/116; % leaves 6 points for the walls at 128^2
    params.Ly   = 128*2.0/116; % leaves 6 points for the walls at 128^2
    params.h    = (params.Lx - 2.0) / 2.0;
    params.nu   = 1/1000;
    params.T_end = 0.40;
    params.x     = params.Ly*(0:params.nx-1)/(params.nx) - params.h;
    params.y     = params.Lx*(0:params.ny-1)/(params.ny) - params.h;
    
    
    % wavenumbers + grid
    geometry_wavenumbers();
end
