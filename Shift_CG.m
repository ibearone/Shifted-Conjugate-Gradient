clc;clear;close all;
%%% Setting Parameters%%%
dim=500 %%%dimention of matrix
%%%%%
 A= diag(randi(1,dim,1))
 R = sprandsym(dim,0.01)
 A=A+R

B=randi(10,dim,1)


nmax=4
sigma=[2 1 0.5 0.1 0]
for n=1:nmax
[x{n},norm_r{n},norm_p{n}]=Shift_ConG(A,B,sigma(n),0.1^5)
end
for n=1:nmax
subplot(1,2,1)
plot(norm_r{n},'LineWidth',1)
hold on;
end
legend('\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')
xlabel Steps
ylabel |r|

set(gca, 'YScale', 'log')
for n=1:nmax
    
subplot(1,2,2)
plot(norm_p{n},'LineWidth',1)
hold on;
end
legend('\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')
xlabel Steps
ylabel |p|

set(gca, 'YScale', 'log')
figure;
for n=1:nmax
    X{n}=inv(A+sigma(n)*eye(dim))*B
for k=1:length(x{n})
def(k,n)=norm(x{n}{k}-X{n})
end
plot(def(:,n),'LineWidth',1)
hold on;
end
legend('\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')
legend('|r|','|x-X|')
xlabel Steps
ylabel |x-X|

set(gca, 'YScale', 'log')