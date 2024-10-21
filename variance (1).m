clear all;
clc;
m=4;          % number of clusters
rho=[0.02, 0.0634, 0.0765, 0.1877];
rho1= rho/2; % ICC

sigma= [0.025,0.042,0.044,0.047];
% rho=[0.06, 0.06, 0.06, 0.06];
% rho1= [0,0,0,0]; % ICC
% sigma= [1,1,1,1];

d=200
T_star=368-d;

c=[1.5, 2, 2.25,2.25];
%c=[2, 2, 2, 2];

[n_opt_exact, n_bal_exact, var_unequal, var_bal] = var_cost_design(m,rho,rho1,sigma,c,T_star);
n_opt_exact
n_bal_exact
%n_bal_1_exact
eff =  var_bal/var_unequal 
 
 N_unequal=sum(n_opt_exact)
 N_bal=m*n_bal_exact


 function[n_opt_exact, n_bal_exact, var_unequal, var_bal] = var_cost_design(m,rho,rho1,sigma,c,T_star)

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
n_opt_exact;


%optimal design
n_bal = T_star/sum(c);
n_bal_exact=round(n_bal);

%variance calculation


for j=1:1:m
   V(j) = (n_opt_exact(j))/((sigma(j)^2)*(1+((n_opt_exact(j)-1)*(rho(j))-rho1(j))));
end

for j=1:1:m
V_bal(j) =  (n_bal_exact)/((sigma(j)^2)*(1+((n_bal_exact-1)*(rho(j))-rho1(j))));
end

var_unequal= 2/sum(V);
var_bal = 2/sum(V_bal);

end

    
