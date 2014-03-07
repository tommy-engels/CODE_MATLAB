function [err,u_new]=time_step(eps,method,dt_fixed,nx)
    % for rk2, skip this test if it would be unstable
    if (strcmp(func2str(method),'rk2'))
        if (dt_fixed > eps)
            err = 1;
            u_new = zeros(nx+6,nx+6);
            return
        end
    end
        

    addpath(genpath('./lib_spectral_matlab/'))
    addpath(genpath('./lib_finite_differences_matlab/'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    dx = 2*pi / (nx-1);
    % grid in the "fluid"
    params.xf = 2*pi*(0:nx-1)/(nx-1);
    params.x = [-3*dx -2*dx -dx params.xf 2*pi+dx 2*pi+2*dx 2*pi+3*dx];
    params.y = [-3*dx -2*dx -dx params.xf 2*pi+dx 2*pi+2*dx 2*pi+3*dx];
    params.nx           = length(params.x);
    params.ny           = length(params.y);
    params.dx=dx;
    params.dy=dx;
    % the grid is periodic and therefore the last point of x is not the length
    % of the domain, but +dx is.
    params.Lx           = params.x(end)-params.x(1) + dx;
    params.Ly           = params.y(end)-params.y(1) + dx;
    params.eta          = eps;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.nu           = 1/10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.T_end        = 0.5; %was 2.5 for error wrt exact solution
    params.dt           = dt_fixed;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [params.X params.Y] = meshgrid(params.x,params.y);
    params.X            = params.X';
    params.Y            = params.Y';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.mask         = zeros(params.nx,params.ny);
    params.dealias      = zeros(params.nx,params.ny);
    params.us           = zeros(params.nx,params.ny);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create wavenumber matrices (global)
    % params.kx             = (2*pi/params.Lx)*[0:(params.nx/2-1) (-params.nx/2):(-1)]; % Vector of wavenumbers
    % params.ky             = (2*pi/params.Ly)*[0:(params.ny/2-1) (-params.ny/2):(-1)]; % Vector of wavenumbers
    params.kx             = fmodes(params.nx,params.Lx);
    params.ky             = fmodes(params.ny,params.Ly);
    [params.Kx,params.Ky] = meshgrid(params.kx,params.ky);
    params.Kx             = params.Kx';
    params.Ky             = params.Ky';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create the mask
    create_mask();
    % initial condition
    uk_old    = inicond();
    u_old     = cofitxy(uk_old);
    uk_oldold = uk_old*0;
    u_oldold  = u_old*0;

    time = 0;
    it   = 1;

    u_ex = exact_solution();

    % figure
    while (time<params.T_end)
        params.dt = min([ dt_fixed, dt_TIME(time,params.T_end) ]);

        [u_new,uk_new] = method(time,it,u_old,uk_old,u_oldold,uk_oldold);

        u_oldold  = u_old;
        uk_oldold = uk_old;    
        u_old     = u_new;
        uk_old    = uk_new;
        time      = time + params.dt;
        it        = it+1;
        
        % diverged? -> brake
        if (max(max(abs(u_new)))>2) 
            err = 1;
            u_new = zeros(nx+6,nx+6);
            return
        end
        
    end

    err = norm(reshape(u_old-u_ex,[],1))/norm(reshape(u_ex,[],1));

    %    if ( time == params.T_end)
    %        clf
    %        subplot 221
    %        pcolor(X,Y,u1);
    %        farge_color
    %        axis equal
    %        colorbar
    %        shading interp
    %        title(['Temp time=' num2str(time) ' dt=' num2str(params.dt,'%e') ])
    %
    %        subplot 223
    %        pcolor(X,Y,u1-u_ex);
    %        farge_color
    %        axis equal
    %        colorbar
    %        shading interp
    %        title(num2str(err))
    %        drawnow            
    %    end
end