function [e,w]=rls(lambda,M,u,d,delta)
% recursive least squares,rls.
% Call:
% [e,w]=rls(lambda,M,u,d,delta)
%
% Input arguments:
% lambda = constant, (0,1]
% M = filter length, dim 1x1
% u = input signal, dim Nx1
% d = desired signal, dim Nx1
% delta = constant for initializaton, suggest 1e-7.
%
% Output arguments:
% e = estimation error, dim Nx1
% w = final filter coefficients, dim Mx1
% Step1:initialize
% 2017-4-4 14:34:33, Author: Gui
w=zeros(M,1);
P=eye(M)/delta;
u=u(:);
d=d(:);
% input signal length
N=length(u);
% error vector
e=d.';
% Step2: Loop, RLS
for n=M:N
    uvec=u(n:-1:n-M+1);
    e(n)=d(n)-w'*uvec;
    k=lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
    P=lambda^(-1)*P-lambda^(-1)*k*uvec'*P;
    w=w+k*conj(e(n));
end


