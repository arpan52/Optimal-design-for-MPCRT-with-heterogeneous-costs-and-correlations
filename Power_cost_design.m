clear all;
clc
m=4;          % number of clusters
rho=[0.02, 0.0634, 0.0765, 0.1877];
rho1= rho/2; % ICC

sigma= [0.025,0.042,0.044,0.047];
% rho=[0.06, 0.06, 0.06, 0.06];
% rho1= [0,0,0,0]; % ICC
% sigma= [1,1,1,1];
T_star=168;

c=[0.5, 2, 2.5,3];
%c=[2, 2, 2, 2];

[n_opt_exact, var_unequal] = var_cost_design(m,rho,rho1,sigma,c,T_star);

delta = 0.02;
alpha = 0.05; 
N_opt = sum(n_opt_exact)
[P1,P2] = pow_cost_design(del,al,m,var_unequal)



function [P1,P2] = pow_cost_design(delta,alpha,m,var_unequal)
z = norminv(1-(alpha/2));
t1 = tinv(alpha/2,m-1);
t2 = tinv(1-(alpha/2),m-1);
P1 = normcdf(-z+(delta/sqrt(var_unequal))) + normcdf(-z-(delta/sqrt(var_unequal)));
P2 = 1- nctcdf(t2,m-1,delta/sqrt(var_unequal))+nctcdf(t1,m-1,delta/sqrt(var_unequal));
end

function[n_opt_exact, var_unequal] = var_cost_design(m,rho,rho1,sigma,c,T_star)

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
% Exact optimal design 
n_opt_exact=round(n_opt_exact);


%variance calculation


for j=1:1:m
   V(j) = (n_opt_exact(j))/((sigma(j)^2)*(1+((n_opt_exact(j)-1)*(rho(j))-rho1(j))));
end


var_unequal= 2/sum(V);

end
