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
    nlk1 = rhs_up( time, uk_new, u_new, 'no' );
    nlk1 = add_pressure (nlk1);   
    
    % exponential factor:    
    vis = cal_vis_2d( dt );
    uk_tmp = vis .* (uk_new + dt*nlk1 );                
    uk_tmp = dealias_2d ( uk_tmp );
    u_tmp = cofitxy_2d(uk_tmp);
    % END OF EULER STEP
    nlk2 = rhs_up( time+dt, uk_tmp, u_tmp, 'no' );
    nlk2 = add_pressure (nlk2);
    
    % advance in time
    uk_new = vis.*uk_new + 0.5*dt*(nlk1.*vis + nlk2);            
    uk_new = dealias_2d ( uk_new );
    u_new = cofitxy_2d(uk_new);
    %----------------------------------------------------------------------
    % 3rd strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    [u_new, uk_new] = penalization_step(time,0.5*dt,u_new,uk_new);
    pk_new=pk;
    
    % see if we survived
    if (max(max(max(abs(u))))>1e2)
        error('split diverged..')
    end
end




function [u_new,uk_new]=penalization_step(time,dt,u,uk)
    global params    
    eps = params.eta;    
    us = create_us( u );
    
    if (dt>eps)
        dt_local = dt/floor(dt/eps);
    else
        dt_local = dt;
    end
    
    t=time;
    while t<time+dt
        penal(:,:,1) = -params.mask.*(u(:,:,1)-us(:,:,1))/params.eta;
        penal(:,:,2) = -params.mask.*(u(:,:,2)-us(:,:,2))/params.eta;
        
        penalk = coftxy_2d(penal);
        penalk = add_pressure(penalk);
        penalk = dealias_2d(penalk);
        penal = cofitxy_2d(penalk); 
        
        u = u + dt_local*penal;        
        t = t + dt_local;
    end
    
    % final projection    
    uk_new = coftxy_2d(u);
    uk_new = dealias_2d(uk_new);
    u_new = cofitxy_2d(uk_new);
end




% function [u_new,uk_new]=penalization_step(time,dt,u,uk)
%     global params    
%     chi = params.mask;
%     eps = params.eta;    
%     us = create_us( u );
%     u_new = zeros(params.nx,params.ny,2);
%     
%     % compute penalization pressure
%     penal(:,:,1) = -params.mask.*(u(:,:,1)-us(:,:,1))/params.eta;
%     penal(:,:,2) = -params.mask.*(u(:,:,2)-us(:,:,2))/params.eta;
%     pk = poisson ( divergence_2d( coftxy_2d(penal) ) );
%     pk = gradient_2d(pk);
%     p = cofitxy_2d(pk);
%     
%     
%     % solve first subproblem
%     a = dt/eps;    
%     
%     for i=1:2
%         % fluid
%         uf1 = u(:,:,i) - p(:,:,i)*dt;
%         % solid
%         us1 = exp(-chi*a).*(u(:,:,i)-chi.*us(:,:,i)+eps*p(:,:,i)) - chi.*us(:,:,i) +eps*p(:,:,i);
%         % combine both
%         u_new (:,:,i) = chi.*us1 + (1-chi).*uf1;
%     end
%     
%     p=p*0;
%     for i=1:2
%         % fluid
%         uf1 = u(:,:,i) - p(:,:,i)*dt;
%         % solid
%         us1 = exp(-chi*a) .* (u(:,:,i)-chi.*us(:,:,i)-eps*p(:,:,i)) - chi.*us(:,:,i) -eps*p(:,:,i);
%         % combine both
%         u_new2 (:,:,i) = chi.*us1 + (1-chi).*uf1;
%     end
%     
%     fprintf('%e %e\n',max(max(abs(cofitxy(divergence_2d(coftxy_2d(u_new)))))),...
%         max(max(abs(cofitxy(divergence_2d(coftxy_2d(u_new2)))))));
%     
%     % final projection    
%     uk_new = coftxy_2d(u_new);
%     uk_new = dealias_2d(uk_new);
%     uk_new = make_incompressible_2d(uk_new);
%     u_new = cofitxy_2d(uk_new);
% end