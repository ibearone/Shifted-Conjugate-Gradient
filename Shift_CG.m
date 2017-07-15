clc;clear;close all;
%%% Setting Parameters%%%
dim=100 %%%dimention of matrix
%%%%%
A=randi(10,dim)
A=(A+transpose(A))/2

B=randi(10,dim,1)

X=inv(A)*B
nmax=10
for n=1:nmax
    sigma(n)=0.1^n
[x{n},norm_r{n},norm_p{n}]=Shift_ConG(A,B,sigma(n),0.00000000001)

subplot(2,1,1)
plot(norm_r{n},'-o')
hold on;
subplot(2,1,2)
plot(norm_p{n},'-o')
hold on;
end
set(gca, 'YScale', 'log')
figure;
for n=1:nmax
for k=1:length(x{n})
def(k,n)=norm(x{n}{k}-X)
end
plot(def(:,n),'-o')
hold on;
end


set(gca, 'YScale', 'log')