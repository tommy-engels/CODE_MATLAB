%% evaluates the reference solution by lamballais on our grid
% when first calling it, the function saves many things in the global
% structure (like FFTs of data , the actual data, etc) so that this
% must not be done for every grid point. still quite slow.
% the original code was much simpler, and I guess one day we can optimize
% it (just take the list of points "XY" in the original code)

function [u,v,p,vor]=lamballais(x,y)
    global params
    
    XY = [x,y];
    Nxy = 1;

    if (params.lamballais_ready==0)
        load ./case_lamballais/data_lamballais/r.dat
        load ./case_lamballais/data_lamballais/theta.dat
        load ./case_lamballais/data_lamballais/v.dat
        load ./case_lamballais/data_lamballais/w.dat
        load ./case_lamballais/data_lamballais/p.dat
        load ./case_lamballais/data_lamballais/om.dat
        params.lamballais.r = r;
        params.lamballais.theta = theta;
        params.lamballais.v = v;
        params.lamballais.w = w;
        params.lamballais.p = p;
        params.lamballais.om = om;
        %  Fourier transform
        Nt = length(theta);
        params.lamballais.Nt = length(theta);        
        params.lamballais.vf = fft(v,Nt,2);
        params.lamballais.wf = fft(w,Nt,2);
        params.lamballais.pf = fft(p,Nt,2);
        params.lamballais.of = fft(om,Nt,2);   
        
        % Chebyshev transform
        Rinf = r(end); Nr=length(r);
        D0=[];
        vec=(0:1:Nr-1)';
        for j=0:1:Nr-1
            D0 = [D0  cos(j*pi*vec/(Nr-1))];
        end
        params.lamballais.D0i = inv(D0);
                
        % mark we're ready
        params.lamballais_ready = 1;
        fprintf('Lamballais interpolation: data loaded...\n');
    end
    
    % load data from global structure
    r     = params.lamballais.r;
    theta = params.lamballais.theta;
    v     = params.lamballais.v;
    w     = params.lamballais.w;
    p     = params.lamballais.p;
    om    = params.lamballais.om;

    % fourier transform
    Nt = params.lamballais.Nt;
    vf = params.lamballais.vf;
    wf = params.lamballais.wf;
    pf = params.lamballais.pf;
    of = params.lamballais.of;
    
    nn = [(0:Nt/2-1) 0 (-Nt/2+1:-1)];

    % Chebyshev transform
    Rinf = r(end);
    Nr = length(r);
    vec=(0:1:Nr-1)';
    D0i = params.lamballais.D0i;

    for m=1:Nxy
        % polar cordinates
        [TI,RI] = cart2pol(XY(m,1),XY(m,2)); 
        if (RI > Rinf)
             U(m,1)=1;
             V(m,1)=0;
             P(m,1)=0;
             OM(m,1)=0;  % extrapolation with uniform velocity outside
        elseif (RI<.5)
             U(m,1)=0;
             V(m,1)=0;
             P(m,1)=0;
             OM(m,1)=0;  % 0 inside the cylinder 
        else
             % first step: interpolation in azimutal direction
             tmp = repmat(exp(1i*nn*TI),Nr,1);
             v1(:,1) = real(mean( vf(:,:).*tmp, 2));
             w1(:,1) = real(mean( wf(:,:).*tmp, 2));
             p1(:,1) = real(mean( pf(:,:).*tmp, 2));
             o1(:,1) = real(mean( of(:,:).*tmp, 2));

             % second step: interpolation in radial direction
             v2 = cos(vec*acos(1-2*(RI-.5)/(Rinf-.5)))'*( D0i*v1 );
             w2 = cos(vec*acos(1-2*(RI-.5)/(Rinf-.5)))'*( D0i*w1 );
             p2 = cos(vec*acos(1-2*(RI-.5)/(Rinf-.5)))'*( D0i*p1 );
             o2 = cos(vec*acos(1-2*(RI-.5)/(Rinf-.5)))'*( D0i*o1 );

             % velocity in cartesian coordinates
             U(m,1) = v2*cos(TI)-w2*sin(TI);
             V(m,1) = v2*sin(TI)+w2*cos(TI);
             P(m,1) = p2;
             OM(m,1) = o2;
        end
    end

    % return values
    u = U;
    v = V;
    p = P;
    vor = OM;
end