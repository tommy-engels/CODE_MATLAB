function [e1] = time_step (eps,dtt,method)
global params
NAME = 'all.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS_guermond()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('overwriting dt...\n')
params.dt = dtt;
fprintf('overwriting eps...\n')
params.eta = eps;


% create the mask + solid velocity
create_mask();
% initial condition
[u, uk] = inicond();

time = 0; 
it = 1;

pk = poisson( divergence_2d( rhs_up(time,uk,u,'yes')) );

tic
while ( time < params.T_end )
   % determine time step 
   if strcmp(params.dt_smaller_eps,'yes')
       dt = min([ dt_CFL(u), dt_EPS(params.eta), dt_TIME(time,params.T_end) ]);
   else
       dt = min([ dt_CFL(u), dt_TIME(time,params.T_end) ]);
   end
   if strcmp(params.dt_fixed,'yes')
       dt = params.dt;
   end
   
   % scheme
   [u,uk,pk] = method(time,dt,u,uk,pk);   
      
   % iterate..
   time = time + dt;
   it = it+1;
   
   % progress
   if (mod(it,params.iprogress)==0)
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
e1 = error_ref(time, u, vor);

% otherwise you cannot access data
save(NAME)
end



