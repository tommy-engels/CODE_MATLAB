function [u_new, uk_new, pk] = RK2_classic(u,uk,pk,dt)
%% classic 2nd order RK2 scheme with exact integration of the diffusion term
%--------------------------------------------------------------------------
% explicit treatment of the penalization term
% implicit diffusion
% pressure is div(nlk) 
% pk (in/out) is a dummy argument
%--------------------------------------------------------------------------
    global params
    
    if strcmp(params.dt_smaller_eps,'no') && (dt > params.eta)
        warning('RK2_classic: time step restriction violated, result may be unstable');
    end
    
    
    pk = pk*0; % dummy argument
    
    % compute non-linear terms
    nlk = nonlinear (uk,u,'yes');
    nlk = add_pressure (nlk);
    
    % exponential factor:
    vis = cal_vis_2d( dt );
    uk_new = vis .* (uk + dt*nlk );
    uk_new = dealias_2d( uk_new );
    u_new = cofitxy_2d(uk_new);
        
    % end of euler step 
    nlk2 = nonlinear(uk_new, u_new, 'yes');
    nlk2 = add_pressure (nlk2);
    
    % advance in time
    uk_new = vis.*uk + 0.5*dt*(nlk.*vis + nlk2);
    uk_new = dealias_2d( uk_new );
        
    % go to physical space (consistent output)
    u_new = cofitxy_2d(uk_new);
    
    if (max(max(max(abs(u))))>1e2)
        error('traditional diverged..')
    end
end