function uk = inicond()
    global params
    u = sin(params.X).*sin(params.Y).*(1-params.mask);
    uk = fft2( u );
end