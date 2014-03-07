function u_ex = exact_solution()
    global params
    x = params.xf(1:end-1); % attention, the last point is counted twice (periodicity)
    y = params.xf(1:end-1);
    [X,Y]=meshgrid(x,y);
    X = X';
    Y = Y';
    u  = sin(X).*sin(Y);
    uk = fft2(u);
    [nx, ny] = size(X);
    
    kx      = fmodes(nx, 2*pi); % domain length is 28pi but last point in xf isn't
    ky      = fmodes(nx, 2*pi);
    
    [Kx,Ky] = meshgrid(kx,ky);
    Kx      = Kx';
    Ky      = Ky';
    vis = exp(-params.nu*params.T_end*(Kx.^2+Ky.^2) );
    uk_final = vis .* uk;
    uex = cofitxy(uk_final);

    % copy exact solution into bigger domain
    [dummy,i]=min(abs(params.x-0));
    [dummy,j]=min(abs(params.x-2*pi));
    
    u_ex = zeros(params.nx,params.ny);
    u_ex( i:j-1,i:j-1 ) = uex;
end