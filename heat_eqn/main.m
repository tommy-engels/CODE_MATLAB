function main
    clear all
    close all
    
    [err, u]=time_step( 1e-4, @strang_spectral, 1e-3, 64 );
    err
end


% function main 
%     clear all
%     % detailed dt/eta space
%     eps = logspace(-4,-1, 30);
%     dt  = logspace(-4,-1, 30);
%     n = 0;
%     
%     %% attention this is wrong
%     C = 64*sqrt(3e-2)/(2*pi);
%     nx = round( 1 + (2*pi)./(C*sqrt(eps)) );
%     nx = nx + mod(nx,2);
%     
%     methods = {@rk2,@sbdf2,@strang_spectral,@swss_spectral,@lie_spectral,@lie_spectral_reverse,@cn2_integrating,@cn2_fd,@strang_fd};
%     % loop over methods
%     for k = 1:length(methods)
%         names{k} = func2str( methods{k} );        
%         
%         for i=1:length(eps)
%             for j=1:length(dt)
%                 [err(i,j,k), u]=time_step( eps(i),methods{k},dt(j),nx(i) );
%                 n=n+1;
%                 fprintf('%i of %i ; nx=%i\n',n,length(eps)*length(dt)*length(methods),nx(i))
%             end
%         end
%         
%         save heat_strang_const_modified3e-2.mat
%         
%         [EPS,DT]=meshgrid(eps,dt);
%         EPS=EPS';
%         DT=DT';
%         
%         figure
%         contourf(EPS,DT,log10(err(:,:,k)), -4:0.10:-1 )        
%         caxis([-4 -1])
%         set(gca,'yscale','log');    set(gca,'xscale','log');    xlabel('eta');    ylabel('dt');        colorbar;
%         line([1e-4 1e-1],[1e-4 1e-1],'color',[1 1 1])
%         saveas(gca,['heat_' names{k} '_dxproptoeta.eps'],'epsc')
%         title(names{k});
%         
%         save heat_strang_const_modified3e-2.mat
%     end
% end


% function main
%     clear all
%     %% time convergence test
%     dt = logspace(-4,-1,10);
%     nx = 512;
%     eps = 1e-1;
% %     methods = {@rk2,@sbdf2,@strang_spectral,@strang_fd,@cn2_integrating,@lie_spectral, @lie_spectral_reverse,@swss_spectral,@cn2_fd};
% %     methods = {@rk2,@sbdf2,@strang_spectral,@cn2_integrating,@lie_spectral,@lie_spectral_reverse,@swss_spectral};%,@strang_fd,@cn2_fd};
%     methods = {@strang_spectral}
%     % just one reference solution by the explicit scheme..
%     [dummy, u_ref] = time_step( eps,@rk2,min([0.98*eps 1e-4]),nx );
%         
%     % run naming
%     for j = 1:length(methods)
%         names{j} = func2str( methods{j} );
%     end
%     
%     % dt loops
%     for i = 2:length(dt)
%         for j=1:length(methods)
%             [dummy, u] = time_step( eps,methods{j},dt(i),nx );
%             err(i,j) = norm(reshape(u-u_ref,[],1)) / norm(reshape(u_ref,[],1));
%         end
%     end
%     
%     % plot
%     figure
%     loglog(dt, err,'o-')
%     legend( names )
%     title(['eta=' num2str(eps)])    
% end
