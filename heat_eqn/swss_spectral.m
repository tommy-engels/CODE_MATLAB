function [u_new, uk_new] = swss_spectral(time,it,u_old,uk_old,u_oldold,uk_oldold)
    %% use SWSS splitting instead of strang splitting (almost the same result as strang)
    global params
    chi = params.mask;
    us  = params.us;
    %----------------------------------------------------------------------
    % Part A of SWSS: first penalization, then heat
    %----------------------------------------------------------------------
    % first penalization...
    a = params.dt / params.eta;
    u_a  = (u_old-chi.*us).*exp(-chi*a) + chi.*us;
    % ...now heat eqn
    vis = exp(-params.nu*params.dt*(params.Kx.^2+params.Ky.^2) );  
    uk_a = dealias(fft2(u_a));
    uk_a = vis .* uk_a;    
    %----------------------------------------------------------------------
    % Part B of SWSS: first heat, then penalization
    %----------------------------------------------------------------------
    % first heat...
    vis = exp(-params.nu*params.dt*(params.Kx.^2+params.Ky.^2) );  
    uk_b = dealias(uk_old);
    uk_b = vis .* uk_b;
    u_b = cofitxy( uk_b );    
    % ...then penalization
    a = params.dt / params.eta;
    u_b  = (u_b-chi.*us).*exp(-chi*a) + chi.*us;
    uk_b = dealias( fft2(u_b) );
    %----------------------------------------------------------------------
    % Combine steps A and B
    %----------------------------------------------------------------------
    uk_new = 0.5 * dealias(uk_a + uk_b);
    u_new = cofitxy( uk_new );
end