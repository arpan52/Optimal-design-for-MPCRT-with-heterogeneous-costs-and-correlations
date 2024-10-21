clear all;
clc;

T_star=668;
m=4;
delta = 0.02;
al = 0.05;
t1 = tinv(al/2,m-1);
t2 = tinv(1-(al/2),m-1);

rho=[0.02, 0.0634, 0.0765, 0.1877];
rho1= rho/2; % ICC

sigma= [0.025,0.042,0.044,0.047];



c=[1.5, 2, 2.25,2.25];

[n_opt_exact, n_bal_exact] = var_cost_design(m,rho,rho1,sigma,c,T_star);  % Locally optimal design

% Uniform Bayesian
a1=[0,0,0,0];
b1=[rho(1)+0.2,rho(2)+0.2,rho(3)+0.2,rho(4)+0.2];
a2=[0,0,0,0];
b2= [rho1(1)+0.2,rho1(2)+0.2,rho1(3)+0.2,rho1(4)+0.2];
x0=[T_star/sum(c),T_star/sum(c),T_star/sum(c),T_star/sum(c)];
A=[];
b=[];
Aeq=[c(1),c(2),c(3),c(4)];
beq=T_star;
lb=[1,1,1,1];
ub=[500,500,500,500];

[n_uni1] = fmincon(@(n)var_design(n,a1,b1,a2,b2,m,sigma),x0,A,b,Aeq,beq,lb,ub);
n_uni = round(n_uni1); % Uniform Bayesian Design

% beta Bayesian
y1=[4,4,10,6];
z1=[90,90,70,20];
y2=[4,4,10,6];
z2= [90,90,70,20];
 x0=[40,30,20,10];

% y1=[4,4,4,4];
% z1=[8,8,8,8,8];
% y2= [4,4,4,4];
% z2=[8,8,8,8,8];

%eff_beta = var_beta_design([21,21,21,21],m,sigma,a1,a2,b1,b2)
[n_beta1] = fmincon(@(n)var_beta_design(n,m,sigma,y1,z1,y2,z2),x0,A,b,Aeq,beq,lb,ub);
n_beta = round(n_beta1);  % Beta Bayesian Design


n_minmax = round([84.1747   88.8480   86.0862   75.7102]);

for i=1:1:10e3
r(1) = unifrnd(0,rho(1)+0.2);
r(2) = unifrnd(0,rho(2)+0.2);
r(3) = unifrnd(0,rho(3)+0.2);
r(4) = unifrnd(0,rho(4)+0.2);
r1(1) = unifrnd(0,(rho(1)/2)+0.2);
r1(2) = unifrnd(0,(rho(2)/2)+0.2);
r1(3) = unifrnd(0,(rho(3)/2)+0.2);
r1(4) = unifrnd(0,(rho(4)/2)+0.2);

for j=1:1:m
   V(j) = (n_opt_exact(j))/((sigma(j)^2)*(1+((n_opt_exact(j)-1)*(r(j))-r1(j))));
end

for j=1:1:m
V_uni(j) =  (n_uni)/((sigma(j)^2)*(1+((n_uni-1)*(r(j))-r1(j))));
end

for j=1:1:m
V_beta(j) =  (n_beta)/((sigma(j)^2)*(1+((n_beta-1)*(r(j))-r1(j))));
end

for j=1:1:m
V_minmax(j) =  (n_minmax)/((sigma(j)^2)*(1+((n_minmax-1)*(r(j))-r1(j))));
end

var_unequal= 2/sum(V);
var_uni = 2/sum(V_uni);
var_beta = 2/sum(V_beta);
var_minmax = 2/sum(V_minmax);

%n_bal_1_exact
pwr_opt(i) = 1- nctcdf(t2,m-1,delta/sqrt(var_unequal))+nctcdf(t1,m-1,delta/sqrt(var_unequal));
pwr_uni(i) = 1- nctcdf(t2,m-1,delta/sqrt(var_uni))+nctcdf(t1,m-1,delta/sqrt(var_uni));
pwr_beta(i) = 1- nctcdf(t2,m-1,delta/sqrt(var_beta))+nctcdf(t1,m-1,delta/sqrt(var_beta));
pwr_minmax(i) = 1- nctcdf(t2,m-1,delta/sqrt(var_minmax))+nctcdf(t1,m-1,delta/sqrt(var_minmax));

end
hold on
boxplot([pwr_opt', pwr_uni',pwr_beta',pwr_minmax'])
xticklabels({'Local','Uniform','Beta','Min-max'})
ylabel({'Power'})

hold off

 function[n_opt_exact, n_bal_exact] = var_cost_design(m,rho,rho1,sigma,c,T_star)

for i=1:1:m
a(i)=1-rho(i)-rho1(i);
end
for i=1:1:m
x1(i)= (sqrt(c(i)*a(i)))/(sigma(i)*rho(i));
end
x=sum(x1);

for i=1:1:m
y1(i)=(c(i)*a(i))/rho(i);
end
y=sum(y1);
% check for T

for i=1:1:m
  if  T_star<= (x*sigma(i)*sqrt(c(i)*a(i)))-y
  disp("T is not sufficient for")
  disp(i)
  break
  end
end

%optimal design for unequal cost
for i=1:1:m
    n_opt(i)= (T_star-((x*sigma(i)*sqrt(c(i)*a(i)))-y))/(x*sigma(i)*rho(i)*sqrt(c(i)/a(i)));
    n_opt_exact(i)=round(n_opt(i));
end

%optimal design
n_bal = T_star/sum(c);
n_bal_exact=round(n_bal);

%variance calculation

end

    
function[var_uni_prior] = var_design(n,a1,b1,a2,b2,m,sigma)

for j=1:1:m
k(j) = n(j)/(sigma(j)^2*(b1(j)-a1(j))*(b2(j)-a2(j))*(n(j)-1));
k1(j) = 1-a2(j)+(b1(j)*(n(j)-1));
k2(j) = 1-b2(j)+(b1(j)*(n(j)-1));
k3(j) = 1-a2(j)+(a1(j)*(n(j)-1));
k4(j) = 1-b2(j)+(a1(j)*(n(j)-1));


var(j) = k(j)*((k1(j)*(log(k1(j))-1))-(k2(j)*(log(k2(j))-1))-(k3(j)*(log(k3(j))-1))+(k4(j)*(log(k4(j))-1)));
end
 var_uni_prior = -(1/2)*sum(var);
 
 

end



 function[eff_beta] = var_beta_design(n,m,sigma,y1,z1,y2,z2)
for j=1:1:m
    fun = @(r1,r2)((r1.^(y1(j)-1)).*(r2.^(y2(j)-1)).*((1-r1).^(z1(j)-1)).*((1-r2).^(z2(j)-1)))./((sigma(j).^2).*(1+((n(j)-1).*(r1)-r2)));
    Int(j) = integral2(fun,0,1,0,1); 
   eff_beta1(j) = (n(j)./beta(y1(j),z1(j)).*beta(y2(j),z2(j))).*Int(j);
end
eff_beta = -(1/2).*sum(eff_beta1);
 end