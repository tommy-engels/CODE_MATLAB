function [u_new, uk_new, pk_new] = RK2_SPLIT(u,uk,pk,dt)
    % here we try not to project the RHS, but instead use "traditional"
    % approaches to enforce incompressibility at the end of the time step
    global params
    chi = params.mask;
    eps = params.eta; 
    us  = params.us;
      
    %----------------------------------------------------------------------
    % 1st strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    dt1 = 0.5*dt; a = dt1/eps;    
    u_new (:,:,1) = (u(:,:,1)-chi.*us(:,:,1)).*exp(-chi*a) + chi.*us(:,:,1);
    u_new (:,:,2) = (u(:,:,2)-chi.*us(:,:,2)).*exp(-chi*a) + chi.*us(:,:,2);    
    uk_new = dealias_2d ( coftxy_2d ( u_new ) );
        
    %----------------------------------------------------------------------
    % 2nd strang step, RK2 time step for NS
    %----------------------------------------------------------------------
    % compute non-linear terms
    u_new = cofitxy_2d( uk_new );
    nlk1 = nonlinear( uk_new, u_new, 'no' );
    % add old pressure term
    nlk1 = nlk1 - gradient_2d(pk);    
    
    % exponential factor:    
    vis = cal_vis_2d( dt );
    uk_tmp = vis .* (uk_new + dt*nlk1 );                
    uk_tmp = dealias_2d ( uk_tmp );
    u_tmp = cofitxy_2d(uk_tmp);
    % END OF EULER STEP
    nlk2 = nonlinear( uk_tmp, u_tmp, 'no' );
    nlk2 = nlk2 - gradient_2d(pk);
    
    % advance in time
    uk_new = vis.*uk_new + 0.5*dt*(nlk1.*vis + nlk2);            
    uk_new = dealias_2d ( uk_new );
    
    %----------------------------------------------------------------------
    % 3rd strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    dt1 = 0.5*dt; a = dt1/eps;
    u_new = cofitxy_2d( uk_new );    
    u_new (:,:,1) = (u_new(:,:,1)-chi.*us(:,:,1)).*exp(-chi*a) + chi.*us(:,:,1);
    u_new (:,:,2) = (u_new(:,:,2)-chi.*us(:,:,2)).*exp(-chi*a) + chi.*us(:,:,2);
    uk_new = dealias_2d ( coftxy_2d ( u_new ) );
    
    %----------------------------------------------------------------------
    % projection step ( incremental pressure scheme )    
    %----------------------------------------------------------------------
    divuk = divergence_2d( uk_new );    
    qk = poisson( divuk );
    uk_new = uk_new - gradient_2d( qk );
    % add pressure increment to old pressure
    pk_new = pk + qk / dt;
    
    %----------------------------------------------------------------------    
    % go to physical space (make sure u and uk output is consistent)
    u_new = cofitxy_2d(uk_new);    
    % see if we survived
    if (max(max(max(abs(u))))>1e2)
        error('split diverged..')
    end
end