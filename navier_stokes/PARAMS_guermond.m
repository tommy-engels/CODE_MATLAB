%% Parameters for guermonds testcase
global params
% geometry of the domain
params.nx = 64;
params.ny = 64;
params.Lx = 2.0;
params.Ly = 2.0;
params.xf = params.Ly*(0:params.nx)/(params.nx)-1.0;
params.yf = params.Lx*(0:params.ny)/(params.ny)-1.0;
params.dx = params.xf(2)-params.xf(1);
params.dy = params.yf(2)-params.yf(1);
N = 8; % add N points left and right
params.x = [ params.dx*(-N:-1)+params.xf(1)  params.xf  params.xf(end)+params.dx*(1:N) ];
params.y = [ params.dy*(-N:-1)+params.xf(1)  params.yf  params.xf(end)+params.dy*(1:N) ];
params.nx = length(params.x);
params.ny = length(params.y);
[params.X,params.Y]=meshgrid_t(params.x,params.y);


% viscosity / reynolds number
params.nu = 1.0;


% initial condition
params.inicond = 'quiescent';


% time parameters
params.T_end = 0.5;
params.CFL = 0.1;
params.iplot = 100;
params.iprogress = 500;
params.dt_smaller_eps='yes';


% penalization
params.ipenalization = 'yes';
params.imoving = 'no';
params.imask = 'guermond';
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

% dealiasing
params.dealias = zeros(params.nx,params.ny);