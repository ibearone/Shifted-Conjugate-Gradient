clc;clear;close all;
%%% Setting Parameters%%%
dim=1000 %%%dimention of matrix
%%%%%
 A= diag(randi(10,dim,1))
 R = sprandsym(dim,0.001)
 A=A+R

B=randi(1,dim,1)


nmax=6
sigma=[10 2 1 0.5 0.1 0]
for n=1:nmax
[x{n},norm_r{n},norm_p{n}]=Shift_ConG(A,B,sigma(n),0.1^5)
end
for n=1:nmax
[x_CG{n},norm_r_CG{n}]=CG_LE(A+sigma(n)*eye(dim),B,0.1^5)
end
%%%|r|
for n=1:nmax
subplot(2,2,1)
plot(norm_r{n},'LineWidth',1)
hold on;
end
legend('\sigma=10','\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')
xlabel Steps
ylabel |r|

set(gca, 'YScale', 'log')
for n=1:nmax
subplot(2,2,2)
plot(norm_r_CG{n},'o--','LineWidth',1)
hold on;
end
legend('\sigma=10','\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')
xlabel Steps
ylabel |r|

set(gca, 'YScale', 'log')
%%%%
% figure;
% for n=1:nmax
%     
% subplot(1,3,2)
% plot(norm_p{n},'LineWidth',1)
% hold on;
% end
% legend('\sigma=10','\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')
% xlabel Steps
% ylabel |p|
% 
% set(gca, 'YScale', 'log')


%%%|x-X|
subplot(2,2,3)
for n=1:nmax
    X{n}=inv(A+sigma(n)*eye(dim))*B
for k=1:length(x{n})
def(k,n)=norm(x{n}{k}-X{n})

end
for k=1:length(x_CG{n})
def_CG(k,n)=norm(x_CG{n}{k}-X{n})
end
plot(def(:,n),'LineWidth',1)
hold on;

end
legend('\sigma=10','\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')

xlabel Steps
ylabel |x-X|

set(gca, 'YScale', 'log')
subplot(2,2,4)
for n=1:nmax
plot(def_CG(:,n),'o--','LineWidth',1)
hold on;
end
legend('\sigma=10','\sigma=2','\sigma=1','\sigma=0.5','\sigma=0.1','\sigma=0','location','best')

xlabel Steps
ylabel |x-X|

set(gca, 'YScale', 'log')


figure;
for n=1:nmax
    subplot(2,3,n)
spy(A+sigma(n)*eye(dim))
end