function e = error_ref(vor, u)
    global params
    
    load('Re200.mat')
    vort_romain=double(v);
    [X_ler,Y_ler]=meshgrid(x,y);
    X_ler = X_ler';
    Y_ler = Y_ler';

    
    % cavity gemeometry:
    vort_interp = interp2 ( params.X', params.Y', vor', X_ler, Y_ler, 'spline' );    
   
    figure
    v=10:10:100;
    contour(X_ler,Y_ler,vort_romain,v, 'color',[0 0 0]);
    hold on
    contour(X_ler,Y_ler,vort_interp,v, 'color',[1 0 0]);
    hold on
    contour(params.X,params.Y,params.mask)
    
    
    e = (norm( reshape(vort_interp-vort_romain,[],1) ) ) / (norm(reshape(vort_romain,[],1)));
end
