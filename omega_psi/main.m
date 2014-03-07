function main
close all
clear all
clc
    addpath(genpath('../common_spectral/'))
    addpath(genpath('../common_active/'))
    addpath(genpath('../common_misc/'))
    addpath(genpath('../common_finite_differences/'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global params
params.CASE = 'guermond';
params.nx           = 33;
params.ny           = 33;
params.Lx   = 2;
params.Ly   = 2;
params.eta          = 1.0e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.nu           = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.dt_fixed     = 0; % fix the time step or adapt it dynamically
params.CFL          = 0.05;
params.T_end        = 0.4;
params.dt           = 1e-4;%min(1e-1,params.eta);
params.iplot        = 150; % plot every iplot time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid
% % % params.x = params.Ly*(0:params.nx-1)/(params.nx)-1.0;
% % % params.y = params.Lx*(0:params.ny-1)/(params.ny)-1.0;
% % % params.dx           = params.x(2)-params.x(1);
% % % params.dy           = params.y(2)-params.y(1);

params.xf = params.Ly*(0:params.nx)/(params.nx)-1.0;
params.yf = params.Lx*(0:params.ny)/(params.ny)-1.0;

params.dx           = params.xf(2)-params.xf(1);
params.dy           = params.yf(2)-params.yf(1);
% add this many points to omega_f
N = 8;
params.x = [ params.dx*(-N:-1)+params.xf(1)...
    params.xf...
    params.xf(end)+params.dx*(1:N) ];

params.y = [ params.dy*(-N:-1)+params.xf(1)...
    params.yf...
    params.xf(end)+params.dy*(1:N) ];
params.nx = length(params.x);
params.ny = length(params.y);
    params.Lx           = params.x(end)-params.x(1) + params.dx;
    params.Ly           = params.y(end)-params.y(1) + params.dx;


[params.X params.Y] = meshgrid(params.x,params.y);
params.X            = params.X';
params.Y            = params.Y';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.mask         = zeros(params.nx,params.ny);
params.dealias      = zeros(params.nx,params.ny);
params.us           = zeros(params.nx,params.ny,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create wavenumber matrices (global)
params.kx             = (2*pi/params.Lx)*[0:(params.nx/2-1) (-params.nx/2):(-1)]; % Vector of wavenumbers
params.ky             = (2*pi/params.Ly)*[0:(params.ny/2-1) (-params.ny/2):(-1)]; % Vector of wavenumbers
[params.Kx,params.Ky] = meshgrid(params.kx,params.ky);
params.Kx             = params.Kx';
params.Ky             = params.Ky'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = params.X;
Y = params.Y;

% initial condition
vork_old = Inicond();

% create the mask
create_mask();

time = 0;
it   = 1;
figure

while (time<params.T_end)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    if (time < 0.20)
%        [vork_new, params.dt] = rk2(time,vork_old);
%    else
       [vork_new, params.dt] = rk2_iter(time,vork_old);
%    end
   vork_old = vork_new;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   time = time + params.dt;
   it = it+1;

   if (mod(it,params.iplot)==0)
       clf
       vor = cofitxy(vork_new);
       uk = vor2u(vork_new);
       u=cofitxy_2d(uk);
       pcolor(X,Y,vor);
       colormap(PaletteMarieAll('Vorticity',600,0.3,5,0.25));
       scale = 1.0;
       axis equal
       c = scale*max ( min(min(abs(vor))), max(max(vor)) );
       caxis([-c c])
       colorbar
       shading interp
       hold on
       quiver(X,Y,u(:,:,1),u(:,:,2))
       title(['Vorticity time=' num2str(time) ' dt=' num2str(params.dt,'%e') ])
       drawnow
   end
end 

% e = error_ref( cofitxy(vork_new), vor2u(cofitxy(vork_new)) )
e = error_ref( time, cofitxy_2d(vor2u(vork_new)), 0 )

save all.mat
end























% function e = error_ref(vor, u)
%     global params
%     
%    %% read the data of leriche and generate coordinate axis x,y
%     m_ler = 721;
%     n_ler = 721;
% 
%     fid = fopen('leriche721.vor');
%     Z   = fscanf(fid, '%g %g %g',[3 inf]);
%     fclose(fid);
% 
%     x_ler    = reshape( Z(1,:), m_ler, n_ler ) + 1.0;
%     y_ler    = reshape( Z(2,:), m_ler, n_ler ) + 1.0;
%     vort_ler = reshape( Z(3,:), m_ler, n_ler );
%     
%     x_ler = x_ler(1,:);
%     y_ler = y_ler(:,1);
% 
%     [X_ler,Y_ler] = meshgrid(x_ler,y_ler);
%     X_ler = X_ler';
%     Y_ler = Y_ler';    
%     
%     %% first variant: interpolate my solution to leriche's grid:
%     % interpolate
%     vort_interp = interp2 ( params.X', params.Y', vor, X_ler, Y_ler, 'spline' );    
%     e_c2l = sqrt( trapz(y_ler, trapz(x_ler,(vort_interp-vort_ler).^2))  ) / sqrt( trapz( y_ler,trapz(x_ler,vort_ler.^2)) );
%     
%     %% second variant: interpolate leriche's solution to the cartesian grid
%     mask=params.mask;
%     vor_ler_interp = interp2 ( X_ler', Y_ler', vort_ler, params.X(mask==0), params.Y(mask==0), 'spline' );  
%     e_l2c = norm ( reshape( vor_ler_interp - vor(mask==0),[],1)) /  norm ( reshape( vor_ler_interp,[],1));
%     
%     fprintf('Cart2Leriche= %e  Leriche2Cart= %e\n',e_c2l,e_l2c)
%         
%     figure
%     v = 10:10:300;
%     contour(X_ler,Y_ler,vort_interp,v, 'color',[1 0 0]);
%     hold on
%     contour(X_ler,Y_ler,vort_ler,v, 'color',[0 0 0]);    
%     
%     
%     e = e_c2l;
% end







