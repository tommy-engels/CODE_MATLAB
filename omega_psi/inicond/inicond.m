function vork = Inicond()
    global params
% %     % dipole wall initial condition
% %     x0 = 1.0;
% %     y0 = 1.0;
% %     d  = 0.1; % start position is y0+-d
% %     r0 = 0.1;
% %     we = 299.528385375226;
% %     vor = zeros(params.nx,params.ny);
% %     for ix = 1:params.nx
% %         for iy=1:params.ny
% %             r1 = sqrt( (params.x(ix)-x0)^2 + (params.y(iy)-y0-d)^2 ) / r0;
% %             r2 = sqrt( (params.x(ix)-x0)^2 + (params.y(iy)-y0+d)^2 ) / r0;
% %             vor(ix,iy) = we * (1-(r1)^2) *exp(-(r1)^2) - we * (1-(r2)^2) *exp(-(r2)^2);
% %         end
% %     end
    
    vor = zeros(params.nx,params.ny);
    vork = fft2(vor);
end