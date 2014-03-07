function [u_new, uk_new] = lie_spectral(time,it,u,uk,u_oldold,uk_oldold)
    %% first order lie-splitting with spectral discretization
    global params
    % first the penalty term
    chi = params.mask; us=params.us;
    dt1 = params.dt; 
    a = dt1/params.eta;    
    u_new  = (u-chi.*us).*exp(-chi*a) + chi.*us;
    uk_new = fft2( u_new );
    uk_new = dealias( uk_new );        
    % then the laplacian
    vis = exp(-params.nu*params.dt*(params.Kx.^2+params.Ky.^2) );   
    uk_new = vis .* uk_new;
    uk_new = dealias( uk_new );
    u_new = cofitxy( uk_new );    
end