function D=D22p(n,h)
%% second order second derivative periodic
D=   diag(ones(n,1)*(-2)) ...
    +diag(ones(n-1,1), 1) ...
    +diag(ones(n-1,1),-1);

D(1,n)=1;


D(n,1)=1;

D=sparse(D/(h^2));

end
