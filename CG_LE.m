%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Conjugate Gradient Method Solving Linear Equation
%%% Input: A*x=B to solve x 
%%%        tolerance norm of the residual value
%%%
%%% Output: x
%%% norm_r: norm of residual
%%%
%%% Guanxiong Qu
%%% 2017/07/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,norm_r]=CG_LE(A,B,tolerance)
%%%Initialization
dim=length(A) %%%dimention of matrix

%%% preparation
%x{1}=randi(10,dim,1)
%x{1}=x{1}/norm(x{1})
x{1}=zeros(dim,1)
p{1}=B-A*x{1}
r{1}=p{1}

%%% Initial condition for loop
k=1
loops=k+1
while k<loops
    a(k)=r{k}'*r{k}/(p{k}'*A*p{k})
    x{k+1}=x{k}+a(k)*p{k}
    r{k+1}=r{k}-a(k)*A*p{k}
    b(k)=r{k+1}'*r{k+1}/(r{k}'*r{k})
    p{k+1}=r{k+1}+b(k)*p{k}
    
    norm_r(k)=norm(r{k})
    %%%% Checking convergence
    if norm_r(k)>tolerance
        loops=loops+1
    else
        loops=loops
    end
    k=k+1
end

end 
