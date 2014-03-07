function [u_new, uk_new] = strang_spectral(time,it,u,uk,u_oldold,uk_oldold)
    %% strang splitting for heat eqn, all terms exact, the time error is just the splitting error
    global params
    chi = params.mask; us=params.us;
    %----------------------------------------------------------------------
    % 1st strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    dt1 = 0.5*params.dt; 
    a = dt1/params.eta;    
    u_new  = (u-chi.*us).*exp(-chi*a) + chi.*us;
    uk_new = fft2( u_new );
    uk_new = dealias( uk_new );        
    %----------------------------------------------------------------------
    % 2nd strang step, solve heat eqn for an entire time step
    %----------------------------------------------------------------------
    % exponential factor:    
    vis = exp(-params.nu*params.dt*(params.Kx.^2+params.Ky.^2));    
    uk_new = vis .* uk_new;
    uk_new = dealias( uk_new );
    u_new = cofitxy( uk_new );    
    %----------------------------------------------------------------------
    % 3rd strang step, half a time step for penalization term
    %----------------------------------------------------------------------
    dt1 = 0.5*params.dt; 
    a = dt1/params.eta;  
    u_new  = (u_new-chi.*us).*exp(-chi*a) + chi.*us;
    uk_new = dealias ( fft2(u_new) );
    u_new  = cofitxy( uk_new );
end