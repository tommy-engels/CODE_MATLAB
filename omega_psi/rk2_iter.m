function [vork_new,dt] = rk2_iter(time, vork)
    global params  
    %% time stepping WITHOUT penalization
    % compute non-linear terms
    [nlk,dt] = rhs_omega(time,vork,'no');
    % advance in time
    vis = exp(-params.nu*dt*(params.Kx.^2+params.Ky.^2) );    
    vork_new = vis.*(vork + dt*nlk );                
    % Dealiasing    
    vork_new = dealias(vork_new);
    % end of euler step    
    [nlk2,dt_dummy] = rhs_omega(time+dt,vork_new,'no');     
    % advance in time
    vork_new = vis.*vork + 0.5*dt*(nlk.*vis + nlk2);        
  
    %% iteration
    % fetch velocity
    uk = vor2u(vork_new);
    uk = mean_flow_forcing( uk );
    uf = cofitxy_2d(uk);
    
    % now u is the velocity field after a full time step EE1 w/o
    % penalization
    
    % velocity "error"
    v0 = params.us-uf;
    gamma   = zeros(params.nx,params.ny);
    u_gamma = zeros(params.nx,params.ny,2);
    
    dxdy = params.dx*params.dy;
    err_old =  dxdy*sum(sum( ...
           abs(params.mask.*(v0(:,:,1)-u_gamma(:,:,1)))...
          +abs(params.mask.*(v0(:,:,2)-u_gamma(:,:,2)))...
          ));
    err_new = 0;
    
    it = 0;
    delta = Inf;
    
%     while (abs(err) > 5.0e-4) && (it<1000)
    while ( delta > 1.01) && (it<1000) 
%     for it=1:50


       uu(:,:,1) = fft2(params.mask.*(v0(:,:,1) - u_gamma(:,:,1)));
       uu(:,:,2) = fft2(params.mask.*(v0(:,:,2) - u_gamma(:,:,2)));
       
       gamma = gamma + vorticity_2d( uu );

       
       u_gamma = vor2u ( gamma );
%        u_gamma = mean_flow_forcing( u_gamma );
       u_gamma  = cofitxy_2d( u_gamma );
       
       err_new = dxdy*sum(sum( ...
           +abs(params.mask.*(v0(:,:,1)-u_gamma(:,:,1)))...
           +abs(params.mask.*(v0(:,:,2)-u_gamma(:,:,2)))...
           ));
       delta = (err_old/err_new);
%        fprintf('new=%e old=%e r=%e\n',err_new,err_old,delta)
       err_old=err_new;
       it=it+1;
    end
%     fprintf('---%i\n',it)
    fprintf('time=%e dt=%e iters=%i terminerr=%e\n',time, dt,it,err_old)
    vork_new = vork_new + gamma;
end