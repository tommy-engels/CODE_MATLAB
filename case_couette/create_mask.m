function create_mask
    global params

    X = params.X;
    Y = params.Y;
    x0 = params.Lx/2;
    y0 = params.Ly/2;
    
    % allocate
    params.mask = zeros(params.nx,params.ny);
    params.us   = zeros(params.nx,params.ny,2);
    
    %-----------------------------------
    % COUETTE FLOW (TWO COAXIAL CYLINDERS)
    %-----------------------------------
    % define phi function
    phi1 = (sqrt((X-x0).^2 + (Y-y0).^2) - params.R1) ;
    phi2 =-(sqrt((X-x0).^2 + (Y-y0).^2) - params.R2) ;
    phi = min(phi1,phi2);
    params.mask ( phi <= 0 ) = 1;
    
    % solid velocity
    R = sqrt( (params.X-x0).^2 + (params.Y-y0).^2);
    u = zeros(params.nx,params.ny);
    params.ud_classic = zeros(params.nx,params.ny,2);
    u ( R < 1.2*params.R1 ) = -params.omega*(params.Y( R < 1.2*params.R1 )-y0);
    params.us(:,:,1) = u;
    u ( R < 1.2*params.R1 ) = +params.omega*(params.X( R < 1.2*params.R1 )-x0);
    params.us(:,:,2) = u;
end