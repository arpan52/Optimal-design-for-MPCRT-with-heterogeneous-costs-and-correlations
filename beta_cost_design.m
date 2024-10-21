clear all;
clc;

T_star=168;
m=4;

rho=[0.02, 0.0634, 0.0765, 0.1877];
rho1= rho/2; % ICC
sigma= [0.025,0.042,0.044,0.047];
c=[1.5, 2, 2.25,2.25];

% beta Bayesian
y1=[4,4,10,6];
z1=[90,90,70,20];
y2=[4,4,10,6];
z2= [90,90,70,20];
x0=[40,30,20,10];
A=[];
b=[];
Aeq=[c(1),c(2),c(3),c(4)];
beq=T_star;
lb=[1,1,1,1];
ub=[500,500,500,500];

%eff_beta = var_beta_design([21,21,21,21],m,sigma,a1,a2,b1,b2)
[n_beta1] = fmincon(@(n)var_beta_design(n,m,sigma,y1,z1,y2,z2),x0,A,b,Aeq,beq,lb,ub);
n_beta = round(n_beta1)

 function[eff_beta] = var_beta_design(n,m,sigma,y1,z1,y2,z2)
for j=1:1:m
    fun = @(r1,r2)((r1.^(y1(j)-1)).*(r2.^(y2(j)-1)).*((1-r1).^(z1(j)-1)).*((1-r2).^(z2(j)-1)))./((sigma(j).^2).*(1+((n(j)-1).*(r1)-r2)));
    Int(j) = integral2(fun,0,1,0,1); 
   eff_beta1(j) = (n(j)./beta(y1(j),z1(j)).*beta(y2(j),z2(j))).*Int(j);
end
eff_beta = -(1/2).*sum(eff_beta1);
 end