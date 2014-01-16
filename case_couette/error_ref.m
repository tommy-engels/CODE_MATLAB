function e = error_ref(vor, u)
    global params
        
    uabs    = sqrt(u(:,:,1).^2+u(:,:,2).^2);
    uabs_ex = sqrt(params.u_ex(:,:,1).^2+params.u_ex(:,:,2).^2);
    e= norm(reshape( (1-params.mask).*(uabs-uabs_ex),[],1),2) / ...
        norm(reshape( (1-params.mask).*uabs_ex,[],1),2);
end
