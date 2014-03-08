function time_step
    global params
    %%%%%%%%%%%%%%%%%%%%
    PARAMS_impulse_cylinder()
%     PARAMS_guermond()
%     PARAMS_guermond_periodic()
    %%%%%%%%%%%%%%%%%%%%
    
    method = @rk2_iter;
    
    % initial condition
    vork = inicond();
    
    % create the mask
    create_mask();
    
    time = 0;
    it   = 1;
    figure
    
    while (time<params.T_end)
        [vork, params.dt] = method(time,vork);
        
        time = time + params.dt;
        it = it+1;
        
        if (mod(it,params.iplot)==0)||(time>params.T_end-params.dt)
            clf
            vor = cofitxy(vork);
            uk = vor2u(vork);
            uk = mean_flow_forcing(uk);
            u = cofitxy_2d(uk);
            stream = streamfct( vork );
            pcolor(params.X,params.Y,vor);
            colormap(PaletteMarieAll('Vorticity',600,0.3,5,0.25));
            scale = 1.0;
            axis equal
            c = scale*max ( min(min(abs(vor))), max(max(vor)) );
            caxis([-c c])
            colorbar
            shading interp
            hold on
%             quiver(params.X,params.Y,u(:,:,1),u(:,:,2))
            title(['Vorticity time=' num2str(time) ' dt=' num2str(params.dt,'%e') ])
            hold on
            contour(params.X,params.Y,stream,20,'color','white');
            hold on
            if strcmp(params.imask,'cylinder')
            phi = sqrt((params.X-0.5*params.Lx).^2 + (params.Y-0.5*params.Ly).^2) - 1.0;
            contour(params.X,params.Y,phi,[0 1000],'color','white');
            end
            drawnow
        end
    end
    
    uk = vor2u(vork);
    uk = mean_flow_forcing(uk);
    u = cofitxy_2d(uk);
    vor = cofitxy(vork);
    e = error_ref( time, u, vor );
    
    save(params.name)
end