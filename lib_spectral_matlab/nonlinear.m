function [nlk] = nonlinear (uk, u, penal, diffusion)
    %% [nlk] = nonlinear (uk, u, penal, diffusion)
    % computes the non-linear terms + penalization in Fourier space.    
    % ---------------------------------------------------------------------
    % INPUT:
    % uk: velocity field (2d) in Fourier space
    % u: velocity field (2d) in physical space
    % OPTIONAL:
    % penal [yes]/no: penalization term yes or no
    % diffusion yes/[no]: add explicit diffusion term
    % ---------------------------------------------------------------------
    % OUTPUT:
    % nlk: the non-linear terms in Fourier-space
    % ---------------------------------------------------------------------
    global params
    
    % optional parameters:
    switch nargin
        case 2
            penal='yes';
            diffusion='no';
        case 3
            diffusion='no';
        otherwise
            error('nonlinear: not enough input arguments')
    end
    
    vor = cofitxy( vorticity_2d(uk) );
    
    %% non-linear transport, optionally penalization
    if strcmp(penal,'yes')
        chi = params.mask/params.eta;
        nlk(:,:,1) = fft2( +vor .* u(:,:,2) - chi.*(u(:,:,1)-params.us(:,:,1)) );
        nlk(:,:,2) = fft2( -vor .* u(:,:,1) - chi.*(u(:,:,2)-params.us(:,:,2)) );
    else
        nlk(:,:,1) = fft2( +vor .* u(:,:,2) );
        nlk(:,:,2) = fft2( -vor .* u(:,:,1) );        
    end  
    
    %% add explicit diffusion term, if set
    if strcmp(diffusion,'yes')
        nlk(:,:,1) = nlk(:,:,1) - params.nu *( (params.Kx.^2+params.Ky.^2).*uk(:,:,1) );
        nlk(:,:,2) = nlk(:,:,2) - params.nu *( (params.Kx.^2+params.Ky.^2).*uk(:,:,2) );
    end
    
    %% sponge technology
    if strcmp(params.sponge,'yes')
        vor = vor.*params.masksponge;
        vortk = fft2(vor);
        nlk = nlk - vor2u(vortk);
    end
    
    nlk = dealias_2d ( nlk );
end
