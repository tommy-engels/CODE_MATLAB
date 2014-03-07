function [u_new, uk_new] = cn2_integrating(time,it,u_old,uk_old,u_oldold,uk_oldold)
    global params
    %% CN2 scheme for the laplacian with exact integration of penalization term
    % you were very proud of this; it should be stable and treat the
    % diffusion implicitly while the penalization is treated exactly. Much
    % to your surprise, the error saturates as it does in splitting methods
    chi = params.mask / params.eta;
    dt = params.dt; nu = params.nu;
    % laplace operator
    K2 = -params.Kx.^2 - params.Ky.^2;    
    laplace = cofitxy(nu*K2.*uk_old);
    b = exp(-chi*dt) .* (u_old + laplace*0.5*dt);
    bk = dealias( fft2(b) );
    % again we can exploit the diagonality in Fourier space to "solve" the
    % linear system
    uk_new = bk ./ (1-0.5*dt*K2*nu);    
    u_new = cofitxy(uk_new);
end