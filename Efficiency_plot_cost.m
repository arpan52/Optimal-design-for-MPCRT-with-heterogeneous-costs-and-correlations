clear all;
clc;

for i=1:1:11
    T_star=168+((i-1)*50);

m=4;          % number of clusters
rho=[0.02, 0.0634, 0.0765, 0.1877];
rho1= rho/2; % ICC

sigma= [0.025,0.042,0.044,0.047];



c=[1.5, 2, 2.25,2.25];

[n_opt_exact, n_bal_exact] = var_cost_design(m,rho,rho1,sigma,c,T_star);

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
n_uni = round(n_uni1);



for j=1:1:m
   V(j) = (n_opt_exact(j))/((sigma(j)^2)*(1+((n_opt_exact(j)-1)*(rho(j))-rho1(j))));
end

for j=1:1:m
V_uni(j) =  (n_uni)/((sigma(j)^2)*(1+((n_uni-1)*(rho(j))-rho1(j))));
end

for j=1:1:m
V_bal(j) =  (n_bal_exact)/((sigma(j)^2)*(1+((n_bal_exact-1)*(rho(j))-rho1(j))));
end

var_unequal= 2/sum(V);
var_uni = 2/sum(V_uni);
var_bal = 2/sum(V_bal);

%n_bal_1_exact
eff_opt(i) =  var_bal/var_unequal;
eff_uni(i) =  var_bal/var_uni; 




end


x=168:50:668;
hold on
plot(x,eff_opt,'o-')
xlabel('$T$','Interpreter','Latex') 
ylabel('$Eff(n,n_B)$','Interpreter','Latex') 
plot(x,eff_uni,"^-")
legend('$Eff(\tilde{n},n_B)$','$Eff(n_U,n_B)$','Interpreter','Latex')
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