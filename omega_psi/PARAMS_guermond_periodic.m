%% Parameters for guermonds testcase
global params
% geometry of the domain
params.nx = 64;
params.ny = 64;
params.Lx = 2.0;
params.Ly = 2.0;
params.x = params.Ly*(0:params.nx-1)/(params.nx)-1.0;
params.y = params.Lx*(0:params.ny-1)/(params.ny)-1.0;
params.dx = params.x(2)-params.x(1);
params.dy = params.y(2)-params.y(1);
[params.X,params.Y]=meshgrid_t(params.x,params.y);


% viscosity / reynolds number
params.nu = 1.0;

% initial condition
params.inicond = 'quiescent';

% error computation
params.error = 'guermond';

% time parameters
params.T_end = 0.5;
params.CFL = 0.1;
params.iplot = 100;
params.iprogress = 500;
params.dt_smaller_eps='yes';
params.dt_fixed = 'yes';
params.dt = 1e-4;


% penalization
params.ipenalization = 'no';
params.imoving = 'no';
params.imask = 'empty';
params.eta = 1e-3;
params.active = 'passive';


% sponge
params.sponge='no';

% forcing
params.forcing = 'guermond';

% wavenumbers
params.kx = fmodes( params.nx, params.Lx ); 
params.ky = fmodes( params.ny, params.Ly );
[params.Kx,params.Ky] = meshgrid_t(params.kx,params.ky);

% mean flow
params.meanflow = 'none';

% dealiasing
params.dealias = zeros(params.nx,params.ny);