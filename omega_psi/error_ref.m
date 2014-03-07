function e = error_ref(time, u, vor)
global params
e = 0;

switch params.error
    case 'guermond'
        u_ex = guermond_ex(time);        
        uabs_ex = (1-params.mask).*sqrt( u_ex(:,:,1).^2 + u_ex(:,:,2).^2 );
        uabs    = (1-params.mask).*sqrt( u(:,:,1).^2 + u(:,:,2).^2 );

        e = norm(reshape(uabs-uabs_ex,[],1),2)/norm(reshape(uabs_ex,[],1),2);
        
        figure
        pcolor(params.X,params.Y,abs(uabs-uabs_ex))
        colorbar
        colormap('gray')
        shading interp        
        
        figure
        quiver(params.X,params.Y,u(:,:,1),u(:,:,2))
        hold on
        quiver(params.X,params.Y,u_ex(:,:,1),u_ex(:,:,2))
        colorbar
        shading interp        
        
        figure
        plot (params.x, u(:,5+end/2,1), params.x, u_ex(:,5+end/2,1))
        
    case 'cylinder'
        % interpolate velocity on the cylinder surface
        x0 = 0.5*params.Lx;
        y0 = 0.5*params.Ly;
        ns = 129;
        phi = 2*pi*(0:ns-1)/ns;
        xc = x0 + cos(phi);
        yc = y0 + sin(phi);
        ux = interp2( params.X', params.Y', u(:,:,1)', xc,yc);
        uy = interp2( params.X', params.Y', u(:,:,2)', xc,yc);
        figure 
        subplot 221
        plot(phi, ux)
        subplot 222
        plot(phi,uy)
        subplot 223
        plot(phi,sqrt(uy.^2+ux.^2))
        subplot 224
        pcolor(params.X,params.Y,sqrt(u(:,:,1).^2+u(:,:,2).^2))
        shading interp
        colorbar
        hold on
        plot(xc,yc,'w')
        dphi = phi(2)-phi(1);
        e=sum( sqrt(uy.^2+ux.^2) ) * dphi;
end

fprintf('error=%e\n',e)
end
