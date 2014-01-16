function parameters
global params
params.CASE = 'lamballais';
params.Lx   = 7.50;
params.Ly   = 5.00;
params.R1   = 0.50;
params.R2   = 1.75;
params.R3   = 2.25;
params.nu   = 1/40;
params.T_end = 6;
params.x     = params.Lx*(0:params.nx-1)/(params.nx);
params.y     = params.Ly*(0:params.ny-1)/(params.ny);
params.lamballais_ready=0;


% wavenumbers + grid
    geometry_wavenumbers();
end