function e = error_ref(time, u, vor)
global params
u_ex = guermond_ex(time);

uabs_ex = (1-params.mask).*sqrt( u_ex(:,:,1).^2 + u_ex(:,:,2).^2 );
uabs    = (1-params.mask).*sqrt( u(:,:,1).^2 + u(:,:,2).^2 );

% uabs_ex = u_ex(:,:,1);
% uabs    = u   (:,:,1);

e = norm(reshape(uabs-uabs_ex,[],1),2) ;
f = norm(reshape(uabs_ex,[],1),2);

e=e/f;

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

% save tmp.mat
end
