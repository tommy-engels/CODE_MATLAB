function main
    clear all
%     close all    
    [e] = simulation ( 2e-5, 128, @EE1_adjoint_pressure, 0.3 );
end

function [e1] = simulation ( eps, nx, method, CFL )
global params
NAME = 'all.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;
addpath(genpath('./lib_spectral_matlab/'))
addpath(genpath('./lib_active/'))
addpath(genpath('./lib_finite_differences_matlab/'))
addpath(genpath('./case_vortex_pair/'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some parameters are set here
params.nx  = 1.5*nx;
params.ny  = nx;
params.eta = eps;
params.CFL = CFL;
params.dt_smaller_eps='yes';
params.projection='everywhere';
params.counter=1;
params.iplot=100;


% active penalization (only with RK2_dave currently)
params.active='passive'; % 'chantalat', 'dave'
% how to compute the beta field?
params.active_beta='central';% 'upwind', 'spectral'
% time stepper for chantalat's advection-diffusion eqn
params.active_chantalat_stepper='central';


% load general parameters (CASE - specific)
parameters();
% allocate memory
allocate_mem();
% create the mask + solid velocity
create_mask();
% initial condition
[u, uk] = inicond();


time = 0; 
it = 1;
pk = poisson( divergence_2d( nonlinear(uk,u,'yes')) );

tic
while ( time < params.T_end )
   % determine time step 
   if strcmp(params.dt_smaller_eps,'yes')
       dt = min([ dt_CFL(u), dt_EPS(params.eta), dt_TIME(time,params.T_end) ]);
   else
       dt = min([ dt_CFL(u), dt_TIME(time,params.T_end) ]);
   end
   
   % scheme
   [u,uk,pk] = method(u,uk,pk,dt);   
   
   
   % iterate..
   time = time + dt;
   it = it+1;
   % progress
   if (mod(it,500)==0)
       time_left = (params.T_end - time) * toc/time;
       fprintf('%02i%% %s - %s %s nx=%i eps=%2.1e dt=%2.1e CFL=%f\n',round(100*time/params.T_end),...
          secs2hms(time_left),secs2hms(toc),...
          func2str(method),params.nx,params.eta,dt,params.CFL);
   end   
   
   % plot
   if (mod(it,params.iplot)==0)
       clf;
       vor = cofitxy(vorticity_2d(uk));
       pcolor(params.X,params.Y,vor);
       colorbar;
       shading interp; title(num2str(time))
       drawnow
   end   
end 
vor = cofitxy(vorticity_2d(uk));
e1 = error_ref(vor, u)


% otherwise you cannot access data
save(NAME)
end



