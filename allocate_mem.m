function allocate_mem
global params
params.us   = zeros (params.nx,params.ny,2);
params.dealias = zeros(params.nx,params.ny); 
params.mask = zeros (params.nx,params.ny);
params.masksponge = zeros (params.nx,params.ny);
end