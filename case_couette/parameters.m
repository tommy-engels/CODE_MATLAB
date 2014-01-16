function parameters
    global params
    params.CASE = 'couette';
    params.Lx   = 2.5;
    params.Ly   = 2.5;
    params.R1   = 0.5;
    params.R2   = 1.0;
    params.omega= 1.25;
    params.nu   = 1/10;
    params.T_end = 1.0;
    params.x     = params.Ly*(0:params.nx-1)/(params.nx);
    params.y     = params.Lx*(0:params.ny-1)/(params.ny);
    
    % wavenumbers + grid
    geometry_wavenumbers();
    
    
    %% exact solution
    ExactSolution_Couette();
end
