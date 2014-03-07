function create_mask
    global params    
%     params.mask = zeros(params.nx,params.ny);
    params.us = zeros(params.nx,params.ny,2);
    a=-1;
    b=+1;
    [dummy,i]=min(abs(params.x-a));
    [dummy,j]=min(abs(params.y-b));
    
    params.mask = ones(params.nx,params.ny);
    params.mask ( i+1:j-1,i+1:j-1 ) = 0;

    if (sum(sum(params.mask))==0)
        warning('Mask is empty..')
    end
    

% %     [dummy,ixmin] = min(abs(params.x-params.h));
% %     [dummy,ixmax] = min(abs(params.x-params.Ly+params.h));
% %     [dummy,iymin] = min(abs(params.y-params.h));
% %     [dummy,iymax] = min(abs(params.y-params.Lx+params.h));
% %     
% % %     params.mask ( (params.X>=params.Lx-params.h) ) = 1;
% % params.mask ( (params.X<=0.0)|(params.X>=2.0)|(params.Y<=0.0)|(params.Y>=2.0 ) ) = 1;
% %     params.mask ( (params.X<=params.h)|(params.X>=params.Lx-params.h)|(params.Y<=params.h)|(params.Y>=params.Ly-params.h ) ) = 1; 
end