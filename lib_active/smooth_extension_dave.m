function smooth_extension_dave(time, u)
%% create prolongation velocity field, optimized for matlab computation
    function h=h(x)
        if x<1
            h= exp(1-1/(1-x));
        else
            h=0;
        end
    end

    global params
    params.us = zeros(params.nx,params.ny,2);
    
    % delta is boundary layer thickness
    delta = (params.R3-params.R2) /2;
    
    %% compute field of normal derivatives
    ux_x = params.D1x*u(:,:,1) ;
    ux_y = u(:,:,1)*params.D1y' ;
    uy_x = params.D1x*u(:,:,2) ;
    uy_y = u(:,:,2)*params.D1y' ;
    
%     ux_x = cofitxy(1i*params.Kx.*fft2(u(:,:,1)));
%     ux_y = cofitxy(1i*params.Ky.*fft2(u(:,:,1)));
%     uy_x = cofitxy(1i*params.Kx.*fft2(u(:,:,2)));
%     uy_y = cofitxy(1i*params.Ky.*fft2(u(:,:,2)));
    
%     [ux_x,ux_y,uy_x,uy_y] = upwind_differences(u);
    
    % actual beta field (normal derivatives, component-wise)
    beta(:,:,1) = (params.n_x.*ux_x + params.n_y.*ux_y);
    beta(:,:,2) = (params.n_x.*uy_x + params.n_y.*uy_y);
    
    % count points in layers
    npoints=1;
    for ix=1:params.nx
        for iy=1:params.ny
            if ((params.phi(ix,iy) >= -delta) && ( params.phi(ix,iy) <= 0))
                npoints = npoints+1;
            end
        end
    end
    npoints= npoints-1; % counted one too far
    
    % allocate memory
    xi_x = (1:npoints)*0;
    xi_y = (1:npoints)*0;
    iix = (1:npoints)*0;
    iiy = (1:npoints)*0;
    
    
    % store all points in memory (we now know how many, do the same loop)
    ik=1;
    for ix=1:params.nx
        for iy=1:params.ny
            if ((params.phi(ix,iy) >= -delta) && ( params.phi(ix,iy) <= 0))
                xi_x(ik) = params.X(ix,iy) - params.n_x(ix,iy) * params.phi(ix,iy);
                xi_y(ik) = params.Y(ix,iy) - params.n_y(ix,iy) * params.phi(ix,iy);
                iiy(ik) = iy;
                iix(ik) = ix;
                ik=ik+1;
            end
        end
    end
    
    % interpolate all points at once
    [X,Y]=meshgrid(params.x,params.y);
    u_ex_x = interp2(X,Y,params.u_ex(:,:,1)',xi_x,xi_y);
    u_ex_y = interp2(X,Y,params.u_ex(:,:,2)',xi_x,xi_y);
    u_n_x = interp2(X,Y,beta(:,:,1)',xi_x,xi_y);
    u_n_y = interp2(X,Y,beta(:,:,2)',xi_x,xi_y);
    
    
    % loop over points in layers
    for ik=1:npoints
        ix =iix(ik);
        iy= iiy(ik);
        % we are in the extension layers        
        s  = -params.phi(ix,iy) / delta; % s is positive
        
        % daves functions:
%         b0 = 3*h(s)-3*h(2*s)+h(3*s);
%         b1 = (5/2)*h(s)-4*h(2*s)+(3/2)*h(3*s);
        
        % simpler functions:
%         b0 = 2*s^3 -3*s^2 + 1;
        b0 = 3*s^4 - 4*s^3 +1;
        b1 = s^3 -2*s^2 + s;
        
        params.us(ix,iy,1) = u_ex_x(ik)*b0 - delta*b1*u_n_x(ik);
        params.us(ix,iy,2) = u_ex_y(ik)*b0 - delta*b1*u_n_y(ik);
    end
  
    
% %     uboth(:,:,1)=u(:,:,1).*(1-params.mask) + params.mask.*(params.us(:,:,1));
% %     uboth(:,:,2)=u(:,:,2).*(1-params.mask) + params.mask.*(params.us(:,:,2));
% %     
% %     
% %     figure
% %     u1d = uboth(:,end/2,2);
% %     m = params.mask (:,end/2);
% %     plot(params.x,u1d,...
% %          params.x(m==1), u1d(m==1),'o',...
% %          params.x(m~=1) , u1d(m~=1),'x')
% %      
% %        
% %     
% %     figure; 
% % %     pcolor(params.X,params.Y,cofitxy(vorticity_2d(coftxy_2d(uboth))))
% % %     shading interp
% % %     farge_color
% % %     hold on
% %     quiver(params.X,params.Y,uboth(:,:,1), uboth(:,:,2),'color',[0 0 1])
% %     hold on
% %     contour(params.X,params.Y,params.phi,[0 99])
% %     hold on   
% %     quiver(params.X,params.Y,params.us(:,:,1), params.us(:,:,2),'color',[1 0 0])
% % %     hold on   
% % %     quiver(params.X,params.Y,params.u_ex(:,:,1), params.u_ex(:,:,2),'color',[0 1 0])
    
end




function  [ux_x,ux_y,uy_x,uy_y] = upwind_differences(u)
% chantalat uses upwind differences to compute the normal derivative field,
% we shall try to do so as well.
global params
    dx = params.dx;
    dy = params.dy;
    
    ux_x = zeros(params.nx,params.ny);
    ux_y = zeros(params.nx,params.ny);
    uy_x = zeros(params.nx,params.ny);
    uy_y = zeros(params.nx,params.ny);
    
    for ix=2:params.nx-1
        for iy=2:params.ny-1
            % derivative in x-direction
            if (params.phi(ix+1,iy)<params.phi(ix-1,iy))
                ux_x(ix,iy) = (u(ix+1,iy,1)-u(ix,iy,1)) / dx;
                uy_x(ix,iy) = (u(ix+1,iy,2)-u(ix,iy,2)) / dx;
            else
                ux_x(ix,iy) = (u(ix,iy,1)-u(ix-1,iy,1)) / dx;
                uy_x(ix,iy) = (u(ix,iy,2)-u(ix-1,iy,2)) / dx;
            end
            % derivative in y-direction
            if (params.phi(ix,iy+1)<params.phi(ix,iy-1))
                ux_y(ix,iy) = (u(ix,iy+1,1)-u(ix,iy,1)) / dy;
                uy_y(ix,iy) = (u(ix,iy+1,2)-u(ix,iy,2)) / dy;
            else
                ux_y(ix,iy) = (u(ix,iy,1)-u(ix,iy-1,1)) / dy;
                uy_y(ix,iy) = (u(ix,iy,2)-u(ix,iy-1,2)) / dy;
            end
        end
    end
end