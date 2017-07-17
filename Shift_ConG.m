%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Shifted Conjugate Gradient Method Solving Linear Equation
%%% Input: (A+sigma I)*x=B to solve x 
%%%        tolerance norm of the residual value
%%%
%%% Output: x
%%% norm_r: norm of residual
%%%
%%% Guanxiong Qu
%%% 2017/07/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,norm_r,norm_p]=Shift_ConG(A,B,sigma,tolerance)
%%%Initialization
dim=length(A) %%%dimention of matrix
loops=dim*2

%%% preparation
x{1}=zeros(dim,1)
p{1}=zeros(dim,1)
r{1}=randi(1,dim,1)./dim
r{2}=B
a(1)=1
pi(1)=1;pi(2)=1
rho(1)=Inf

%%% Initial condition for loop
k=2
loops=k+1
while k<loops
    rho(k)=r{k}'*r{k}
    b(k-1)=rho(k)/rho(k-1)
    a(k)=rho(k)/(r{k}'*A*r{k}-b(k-1)*rho(k)/a(k-1))
    r{k+1}=(1+a(k)*b(k-1)/a(k-1)).*r{k}-a(k)*A*r{k}-a(k)*b(k-1)/a(k-1)*r{k-1}
    
    
    pi(k+1)=(1+a(k)*sigma)*pi(k)-a(k)*b(k-1)/a(k-1)*(pi(k-1)-pi(k))
    p{k}=1/pi(k).*r{k}+b(k-1)*(pi(k-1)/pi(k))^2.*p{k-1}
    x{k}=x{k-1}+pi(k)/pi(k+1)*a(k).*p{k}
    
    norm_r(k)=norm(r{k})
    norm_p(k)=norm(p{k})
        %%%% Checking convergence
    if norm_r(k)>tolerance & norm_p(k)>tolerance
        loops=loops+1
    else
        loops=loops
    end
    k=k+1
end
end