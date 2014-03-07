function u_ex=guermond_ex(time)
global params
u_ex(:,:,1) = pi*sin(time)*sin(2*pi*params.Y).*(sin(pi*params.X).^2);
u_ex(:,:,2) =-pi*sin(time)*sin(2*pi*params.X).*(sin(pi*params.Y).^2);

u_ex(:,:,1) = u_ex(:,:,1).*(1-params.mask);
u_ex(:,:,2) = u_ex(:,:,2).*(1-params.mask);
end