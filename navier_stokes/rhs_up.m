function [nlk] = rhs_up (time, uk, u, penal, diffusion)
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
        case 3
            penal='yes';
            diffusion='no';
        case 4
            diffusion='no';
    end
    
    nlk = zeros(params.nx,params.ny,2);
    
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
    
    %% forcing (guermond's testcase)
    if strcmp(params.forcing,'guermond')
       X=params.X;
       Y=params.Y;
        
       cost=cos(time);
       sint=sin(time);
       
       Qx = 2*(pi^3)*(sin(time)^2)*cos(pi*X).*(sin(pi*X).^3).*(sin(2*pi*Y).^2)...
             -2*(pi^3)*(sin(time)^2)*(sin(pi*X).^2).*sin(2*pi*X).*(sin(pi*Y).^2).*cos(2*pi*Y);
       Qy = 2*(pi^3)*(sin(time).^2).*(sin(2*pi*X).^2).*cos(pi*Y).*(sin(pi*Y).^3)...
           -2*(pi^3)*(sin(time).^2).*(sin(pi*X).^2).*cos(2*pi*X).*(sin(pi*Y).^2).*sin(2*pi*Y); 
      
       utx = pi*cost*sin(2*pi*Y).*(sin(pi*X).^2);
       uty = pi*cost*(-sin(2*pi*X)).*(sin(pi*Y).^2);
       
       laplx = pi*sint*   (+2*(pi^2)*(2*cos(2*pi*X)-1).*sin(2*pi*Y));
       laply = pi*sint*   (-2*(pi^2)*sin(2*pi*X).*(2*cos(2*pi*Y)-1));
       
       px = sint*  (-pi*sin(pi*X).*sin(pi*Y));
       py = sint*  ( pi*cos(pi*X).*cos(pi*Y));
              
       fx =  utx + px - laplx + Qx; 
       fy =  uty + py - laply + Qy;
       
       fxk = dealias( fft2( fx.*(1-params.mask) ) );
       fyk = dealias( fft2( fy.*(1-params.mask) ) );
       
       nlk(:,:,1) = nlk(:,:,1) + fxk;
       nlk(:,:,2) = nlk(:,:,2) + fyk;
    end
    
    nlk = dealias_2d ( nlk );
end
