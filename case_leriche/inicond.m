function [u0, uk] = inicond()
    global params   
    
    
    x0 = 1.0;
    y0 = 1.0;
    d  = 0.1; % start position is y0+-d
    r0 = 0.1;
    we = 299.528385375226;
    
    vor = zeros(params.nx,params.ny);
    
    for ix = 1:params.nx
        for iy=1:params.ny
            r1 = sqrt( (params.x(ix)-x0)^2 + (params.y(iy)-y0-d)^2 ) / r0;
            r2 = sqrt( (params.x(ix)-x0)^2 + (params.y(iy)-y0+d)^2 ) / r0;
            vor(ix,iy) = we * (1-(r1)^2) *exp(-(r1)^2) - we * (1-(r2)^2) *exp(-(r2)^2);
        end
    end
    stream = poisson( fft2(vor) );
    u0(:,:,1) = cofitxy( -1i*params.Ky.*stream  );
    u0(:,:,2) = cofitxy( +1i*params.Kx.*stream  );
    
    uk = coftxy_2d(u0);
    fprintf('initial field divergence: %12.5e\n',max(max(abs(cofitxy(divergence_2d(uk))))));
end
