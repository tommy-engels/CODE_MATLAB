%% Parameters for guermonds testcase
global params
% geometry of the domain
params.nx = 64;
params.ny = 64;
params.Lx = 2.0;
params.Ly = 2.0;
params.x = params.Lx*(0:params.nx-1)/params.nx;
params.y = params.Lx*(0:params.nx-1)/params.nx;
[params.X,params.Y]=meshgrid_t(params.x,params.y);


% viscosity / reynolds number
params.nu = 1/25;

% initial condition
params.inicond = 'impulsively_x';

% error computation
params.error = 'none';

% time parameters
params.T_end = 0.5;
params.CFL = 0.1;
params.iplot = 100;
params.iprogress = 500;
params.dt_smaller_eps='yes';
params.dt_fixed = 'yes';
params.dt = 1e-4;


% penalization
params.ipenalization = 'yes';
params.imoving = 'no';
params.imask = 'cylinder';
params.eta = 1e-3;
params.active = 'passive';


% sponge
params.sponge='no';

% forcing
params.forcing = 'none';

% wavenumbers
params.kx = fmodes( params.nx, params.Lx ); 
params.ky = fmodes( params.ny, params.Ly );
[params.Kx,params.Ky] = meshgrid_t(params.kx,params.ky);

% dealiasing
params.dealias = zeros(params.nx,params.ny);

% mean flow
params.meanflow = 'impulsively_x';