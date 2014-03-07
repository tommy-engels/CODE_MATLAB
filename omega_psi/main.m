function main
close all
clear all
clc
global params
addpath(genpath('../common_spectral/'))
addpath(genpath('../common_active/'))
addpath(genpath('../common_misc/'))
addpath(genpath('../common_finite_differences/'))
addpath(genpath('./mask/'))
addpath(genpath('./inicond/'))

%%%%%%%%%%%%%%%%%%%%
PARAMS_guermond()
%%%%%%%%%%%%%%%%%%%%

% initial condition
vork_old = inicond();

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







