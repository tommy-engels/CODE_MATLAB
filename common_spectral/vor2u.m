function uk = vor2u(vork)
%% uk = vor2u(vork) computes the velocity field given vorticity in Fourier space
% note irrotational flows, such as an imposed mean flow, are excluded
    global params
    stream = poisson(vork);
    uk(:,:,1) = -1i*params.Ky.*stream;
    uk(:,:,2) = +1i*params.Kx.*stream;
end