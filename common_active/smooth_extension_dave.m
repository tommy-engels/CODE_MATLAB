function u_smooth = smooth_extension_dave( u )
%% create prolongation velocity field, optimized for matlab computation
    function h=h(x)
        if x<1
            h= exp(1-1/(1-x));
        else
            h=0;
        end
    end

    global params
    u_smooth = zeros(params.nx,params.ny,2);
    
    % delta is boundary layer thickness
    delta = (params.R3-params.R2) /2;
    
    % compute field of normal derivatives
    beta = compute_beta( u );
    
    
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
    [X,Y]  = meshgrid(params.x,params.y);
    u_ex_x = interp2(X,Y,params.u_ex(:,:,1)',xi_x,xi_y);
    u_ex_y = interp2(X,Y,params.u_ex(:,:,2)',xi_x,xi_y);
    u_n_x  = interp2(X,Y,beta(:,:,1)',xi_x,xi_y);
    u_n_y  = interp2(X,Y,beta(:,:,2)',xi_x,xi_y);
    
    
    % loop over points in layers
    for ik=1:npoints
        ix = iix(ik);
        iy = iiy(ik);
        % we are in the extension layers        
        s  = -params.phi(ix,iy) / delta; % s is positive
        
        % daves functions:
%         b0 = 3*h(s)-3*h(2*s)+h(3*s);
%         b1 = (5/2)*h(s)-4*h(2*s)+(3/2)*h(3*s);
        
        % simpler functions:
%         b0 = 2*s^3 -3*s^2 + 1;
        b0 = 3*s^4 - 4*s^3 +1;
        b1 = s^3 -2*s^2 + s;
        b0=0;
        u_smooth(ix,iy,1) = u_ex_x(ik)*b0 - delta*b1*u_n_x(ik);
        u_smooth(ix,iy,2) = u_ex_y(ik)*b0 - delta*b1*u_n_y(ik);
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