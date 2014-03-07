function [nlk, dt]=rhs_omega(t, vork, penalization)
    %% computes the non-linear terms + penalization in Fourier space.    
    global params
    uk = vor2u(vork);
    uk = mean_flow_forcing(uk);
    u = cofitxy_2d(uk);
            
    %% determine time step 
    if strcmp(params.dt_fixed,'no')
        dt = min(params.CFL*params.dx/max(max(max(abs(u)))), params.eta);
    else
        dt = params.dt;
    end
    
    %% nonlinear term and penalization
    if strcmp(penalization,'yes')
        % note this is the curl
        penalk(:,:,1) = +1i*params.Kx.*fft2( (params.mask/params.eta).*(u(:,:,2)-params.us(:,:,2)) );
        penalk(:,:,2) = -1i*params.Ky.*fft2( (params.mask/params.eta).*(u(:,:,1)-params.us(:,:,1)) );        
        nlk = fft2( -u(:,:,1).*cofitxy(1i*params.Kx.*vork) - u(:,:,2).*cofitxy(1i*params.Ky.*vork)  )  - (penalk(:,:,1) + penalk(:,:,2));
    else
        nlk = fft2( -u(:,:,1).*cofitxy(1i*params.Kx.*vork) - u(:,:,2).*cofitxy(1i*params.Ky.*vork)  );
    end
    
    %% add source term that forces exact solution in guermonds paper   
    if strcmp(params.forcing,'guermond')
        x = params.X;
        y = params.Y;
        
        F = 4.*pi.^4.*sin(t).^2.*cos(pi.*x).*sin(pi.*x).^3.*cos(2.*pi.*y).*sin(2.*pi.*y)+((-12.*pi.^4.*sin(t)-2.*pi.^2.*cos(t)).*sin(pi.*x).^2+...
            4.*pi.^4.*sin(t).*cos(pi.*x).^2).*cos(2.*pi.*y)-4.*pi.^4.*sin(t).^2.*cos(2.*pi.*x).*sin(2.*pi.*x).*cos(pi.*y).* ...
            sin(pi.*y).^3+(-12.*pi.^4.*sin(t)-2.*pi.^2.*cos(t)).*cos(2.*pi.*x).*sin(pi.*y).^2+4.*pi.^4.*sin(t).*cos(2.*pi.*x).*cos(pi.*y).^2;
        
        nlk = nlk + fft2( (1-params.mask).*F );
    end
    
    %% dealiasing 
    nlk = dealias( nlk );
end
