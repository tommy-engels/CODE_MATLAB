function e = error_ref(vor, u)
    global params
    
   %% read the data of leriche and generate coordinate axis x,y
    m_ler = 721;
    n_ler = 721;

    fid = fopen('leriche721.vor');
    Z   = fscanf(fid, '%g %g %g',[3 inf]);
    fclose(fid);

    x_ler    = reshape( Z(1,:), m_ler, n_ler ) + 1.0;
    y_ler    = reshape( Z(2,:), m_ler, n_ler ) + 1.0;
    vort_ler = reshape( Z(3,:), m_ler, n_ler );
    
    x_ler = x_ler(1,:);
    y_ler = y_ler(:,1);

    [X_ler,Y_ler] = meshgrid(x_ler,y_ler);
    X_ler = X_ler';
    Y_ler = Y_ler';    
    
    %% first variant: interpolate my solution to leriche's grid:
    % interpolate
    vort_interp = interp2 ( params.X', params.Y', vor, X_ler, Y_ler, 'spline' );    
    e_c2l = sqrt( trapz(y_ler, trapz(x_ler,(vort_interp-vort_ler).^2))  ) / sqrt( trapz( y_ler,trapz(x_ler,vort_ler.^2)) );
    
    %% second variant: interpolate leriche's solution to the cartesian grid
    mask=params.mask;
    vor_ler_interp = interp2 ( X_ler', Y_ler', vort_ler, params.X(mask==0), params.Y(mask==0), 'spline' );  
    e_l2c = norm ( reshape( vor_ler_interp - vor(mask==0),[],1)) /  norm ( reshape( vor_ler_interp,[],1));
    
    fprintf('Cart2Leriche= %e  Leriche2Cart= %e\n',e_c2l,e_l2c)
        
    figure
    v = 10:10:300;
    contour(X_ler,Y_ler,vort_interp,v, 'color',[1 0 0]);
    hold on
    contour(X_ler,Y_ler,vort_ler,v, 'color',[0 0 0]);    
    
    
    e = e_c2l;
end
