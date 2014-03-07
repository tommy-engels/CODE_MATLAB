function [u_new, uk_new] = rk2_classic(time,it,u,uk,u_oldold,uk_oldold)
    %% explicit RK2 for penalty term, integrating factor for heat kernel
    global params  
    % compute non-linear terms
    [nlk] = penalty(u,uk);
    % advance in time
    vis = exp(-params.nu*params.dt*(params.Kx.^2+params.Ky.^2) );    
    uk_new = vis.*(uk + params.dt*nlk );                    
    uk_new = dealias(uk_new);
    u_new = cofitxy( uk_new );
    
    % end of euler step    
    [nlk2] = penalty(u_new, uk_new);     
    % advance in time
    uk_new = vis.*uk + 0.5*params.dt*(nlk.*vis + nlk2);        
  
    % Dealiasing
    uk_new = dealias(uk_new);
    u_new = cofitxy( uk_new );
end


function [rhsk] = penalty(u,uk)
    %% penalization term in Fourier space, used only in explicit method
    global params
    rhsk = fft2( -params.mask .* (u-params.us) / params.eta );
    rhsk = dealias( rhsk );        
end
