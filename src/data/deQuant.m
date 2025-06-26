function [xhat,deltaQuant] = deQuant(x,deltaQuant)
% DEQUANT filters a signal according to the minimum quadratic second
% derivative within the quantization step size
%
% [xhat] = deQuant(x,deltaQuant) uses the quantization step size
% deltaQuant to calculate the estimate xhat from the data vector x
%
% [xhat,deltaQuant] = deQuant(x) uses the minimum step size from
% the data vector x as the quantization step size
%
% assumptions
% - constant sample rate
% - constant quantization
% - true signal is smooth

if nargin < 2
    deltaQuant = min(diff(unique(x))); % determine quantization
end

% number of samples
N = length(x);

% abs( xhat_k - x_k ) <= deltaQuant for every sample k:
% xhat - x <= deltaQuant
% -(xhat - x) <= deltaQuant
%-> A*xhat <= b
b = [x+deltaQuant/2 ; -x+deltaQuant/2];
A = [ diag(ones(N,1)) ; diag(-ones(N,1)) ];

% initial guess
xhat0 = x;

% constant objective Hessian
% H =   diag(ones([N-2 1]),2) + diag(ones([N-2 1]),-2)...
%     + diag([-2;-4*ones([N-3 1]);-2],1) + diag([-2;-4*ones([N-3 1]);-2],-1)...
%     + diag([1;5;6*ones([N-4 1]);5;1],0); 
% H=sparse(H);

H=spdiags([ones([N 1]) [-2;-4*ones([N-3 1]);-2;NaN] [1;5;6*ones([N-4 1]);5;1] [NaN;-2;-4*ones([N-3 1]);-2] ones([N 1])],[-2,-1,0,1,2],N,N);

% solver options 
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'HessianFcn',@(x,lambda) H,...
    'Display','iter');

% solve optimization problem
xhat = fmincon(@(xhat) myobjective(xhat,N),xhat0,A,b,[],[],[],[],[],options);

end


function [f,g] = myobjective(xhat,N)

iN = 2:N-1;

f = sum( 0.5*(xhat(iN+1)-2*xhat(iN)+xhat(iN-1)).^2 );

if nargout > 1
    
iNg = 3:N-2;

g = [1*(xhat(1  +2)-2*xhat(1  +1)+xhat(1  ));...
     1*(xhat(2  +2)-2*xhat(2  +1)+xhat(2  ))-2*(xhat(2  +1)-2*xhat(2  )+xhat(2  -1));...
    (1*(xhat(iNg+2)-2*xhat(iNg+1)+xhat(iNg))-2*(xhat(iNg+1)-2*xhat(iNg)+xhat(iNg-1))+1*(xhat(iNg)-2*xhat(iNg-1)+xhat(iNg-2)));...
                                            -2*(xhat(N-1+1)-2*xhat(N-1)+xhat(N-1-1))+1*(xhat(N-1)-2*xhat(N-1-1)+xhat(N-1-2));...
                                                                                    +1*(xhat(N  )-2*xhat(N  -1)+xhat(N  -2))];
% g(isnan(g))=0;

end

end
