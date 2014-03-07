function time_step
    global params
    %%%%%%%%%%%%%%%%%%%%
    PARAMS_impulse_cylinder()
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
               [vork_new, params.dt] = rk2(time,vork_old);
%         %    else
%         [vork_new, params.dt] = rk2_iter(time,vork_old);
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
            pcolor(params.X,params.Y,vor);
            colormap(PaletteMarieAll('Vorticity',600,0.3,5,0.25));
            scale = 1.0;
            axis equal
            c = scale*max ( min(min(abs(vor))), max(max(vor)) );
            caxis([-c c])
            colorbar
            shading interp
            hold on
            quiver(params.X,params.Y,u(:,:,1),u(:,:,2))
            title(['Vorticity time=' num2str(time) ' dt=' num2str(params.dt,'%e') ])
            drawnow
        end
    end
    
    % e = error_ref( cofitxy(vork_new), vor2u(cofitxy(vork_new)) )
    e = error_ref( time, cofitxy_2d(vor2u(vork_new)), 0 )
    
    save all.mat
end