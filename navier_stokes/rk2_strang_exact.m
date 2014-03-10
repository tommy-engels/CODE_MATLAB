%% Runge kutta 2 with strang splitting and separate projection
% the strang splitting advanec the penalizaed navier-stokes eqn without the
% pressure, and the operators
% A(u) = -chi/eta *u
% B(u) = (u.grad)u -nu*laplace(u)
% are evolved serparately. The operator A is solved exactly (an
% exponential)

function [u_new, uk_new, pk_new] = rk2_implicit(time, dt,u,uk,pk)
    global params      
    %----------------------------------------------------------------------
    % 1st strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    [u_new, uk_new] = penalization_step(time,0.5*dt,u,uk);
        
    %----------------------------------------------------------------------
    % 2nd strang step, RK2 time step for NS
    %----------------------------------------------------------------------
    % compute non-linear terms
    u_new = cofitxy_2d( uk_new );
    nlk1 = rhs_up( time, uk_new, u_new, 'no' );
    % add old pressure term
    nlk1 = nlk1 - gradient_2d(pk);    
    
    % exponential factor:    
    vis = cal_vis_2d( dt );
    uk_tmp = vis .* (uk_new + dt*nlk1 );                
    uk_tmp = dealias_2d ( uk_tmp );
    u_tmp = cofitxy_2d(uk_tmp);
    % END OF EULER STEP
    nlk2 = rhs_up( time+dt, uk_tmp, u_tmp, 'no' );
    nlk2 = nlk2 - gradient_2d(pk);
    
    % advance in time
    uk_new = vis.*uk_new + 0.5*dt*(nlk1.*vis + nlk2);            
    uk_new = dealias_2d ( uk_new );
    u_new = cofitxy_2d(uk_new);
    %----------------------------------------------------------------------
    % 3rd strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    [u_new, uk_new] = penalization_step(time,0.5*dt,u_new,uk_new);
    
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


function [u_new,uk_new]=penalization_step(time,dt,u,uk)
    global params    
    chi = params.mask;
    eps = params.eta;     
    
    a = dt/eps;    
    us = create_us( u );
    u_new (:,:,1) = (u(:,:,1)-chi.*us(:,:,1)).*exp(-chi*a) + chi.*us(:,:,1);
    u_new (:,:,2) = (u(:,:,2)-chi.*us(:,:,2)).*exp(-chi*a) + chi.*us(:,:,2);  
    
    uk_new = coftxy_2d(u_new);
    uk_new = dealias_2d(uk_new);
end