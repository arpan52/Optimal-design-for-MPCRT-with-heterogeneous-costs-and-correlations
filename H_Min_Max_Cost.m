

% H-ALGORITHM : 
warning off;

T_star=668;
m=4;
sigma= [0.025,0.042,0.044,0.047];
c=[1.5, 2, 2.25,2.25];


global iteration_num
iteration_num=0;

%Step 0

l=1; B_pi_l = -inf;
h_grid_space = 0.05;
H = 0.0:h_grid_space:1.0;
%theta = [mu;beta;tau;gam];
%theta_vals_l = [ ones(1,1+p+1) ; zeros(1,1+p+1)];
%For prior distribution, assigning equal prob mass to both theta
%prob_vals_l = [0.5, 0.5];
"Initial Prior on Theta : "
theta_vals_l = 0.2*ones(1,8);
prob_vals_l = 1.0;

% Initial Prior PMF is (theta_vals, prob_vals)
"Starting n"; 
n_l = get_optimal_on_the_average_n(theta_vals_l,prob_vals_l,  m,sigma,c,T_star);
B_pi_l = B(n_l, theta_vals_l, prob_vals_l,  m,sigma,c,T_star);

step_1(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,  m,sigma,c,T_star)

function step_1(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,m,sigma,c,T_star)
    global iteration_num;

    theta_argmax_psi = zeros(1,2*m); psi_max = -inf;
    stopping=true;
    
    for i=1:1:1        
        %initialization of theta_nod according to (0.8,1.1) limit
        %theta_nod = 0.3*rand(1,2)+(0.8);
        theta_nod = 0.05*ones(1,2*m);

rho=[0.02, 0.0634, 0.0765, 0.1877];
rho1= rho/2; % ICC
a1=[0,0,0,0,0,0,0,0];
b1=[rho(1)+0.2,rho(2)+0.2,rho(3)+0.2,rho(4)+0.2,rho1(1)+0.2,rho1(2)+0.2,rho1(3)+0.2,rho1(4)+0.2];

options = optimoptions('fmincon','Display','none');
        [this_theta_psi_max, this_psi_max] = fmincon(@(theta) -psi1(n_l,theta, m,sigma,c,T_star), theta_nod, [],[], [], [], a1, b1,[],options);
        this_psi_max = -1.0*this_psi_max;
        if this_psi_max > psi_max
            psi_max = this_psi_max;
            theta_argmax_psi = this_theta_psi_max;
        end
    end
    iteration_num
    iteration_num = iteration_num+1;
    psi_max
    B_pi_l
    if psi_max>B_pi_l
        stopping=false;
    end
              
    if stopping==true
        "Found minmax n"
        n_l
        "Least Favourable Distribution"
        "probs"
        prob_vals_l
        "thetas"
        theta_vals_l
    else
        theta_l = theta_argmax_psi;
        step_2(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,theta_l, m,sigma,c,T_star);
    end
end

function step_2(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,theta_l,m,sigma,c,T_star)
    %delta_l = unit mass to theta_l
    step_3(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,theta_l,  m,sigma,c,T_star);
end

function step_3(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,theta_l,  m,sigma,c,T_star)

n_l1 = n_l;
    largest_B_pi_t_l1 = -inf;
    prob_vals_l1 = prob_vals_l
    theta_vals_l1 = theta_vals_l
    
    theta_vals_l;
    prob_vals_l;
    
    for t_val=1:1:length(H)
        theta_vals_t_l1 = [theta_vals_l; theta_l];
        %"old"
        %prob_vals_l1
        prob_vals_t_l1 = [(1-H(t_val)).*prob_vals_l H(t_val)];
        %"H(t_val)"
        %H(t_val)
        %"new"
        %prob_vals_t_l1
        
        n_t_l1 = get_optimal_on_the_average_n(theta_vals_t_l1, prob_vals_t_l1,  m,sigma,c,T_star);
        B_pi_t_l1 = B(n_t_l1, theta_vals_t_l1, prob_vals_t_l1,  m,sigma,c,T_star);        
        if largest_B_pi_t_l1<B_pi_t_l1
            largest_B_pi_t_l1 = B_pi_t_l1;
            n_l1 = n_t_l1;
            prob_vals_l1 = prob_vals_t_l1;
            theta_vals_l1 = theta_vals_t_l1;
        end
    end
    step_4(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,n_l1, theta_vals_l1, prob_vals_l1, theta_l,  m,sigma,c,T_star)    
end

function step_4(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,n_l1, theta_vals_l1, prob_vals_l1, theta_l,  m,sigma,c,T_star) 

B_pi_l1 = B(n_l1, theta_vals_l1, prob_vals_l1,  m,sigma,c,T_star);
    if B_pi_l1 > B_pi_l
        B_pi_l = B_pi_l1;
        n_l = n_l1;
        "Assigned l1 to l";
        theta_vals_l = theta_vals_l1;
        prob_vals_l = prob_vals_l1; 
        step_1(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,  m,sigma,c,T_star);
    else
        h_grid_space = h_grid_space/2.0;
        H = 0.0:h_grid_space :1.0;
        step_3(H,h_grid_space,theta_vals_l,prob_vals_l,n_l,B_pi_l,theta_l,  m,sigma,c,T_star);
    end
end

% function ret_val = get_least_favourable_argument(n)
% 
%     arg_0 = 10*rand(1,2);options = optimoptions('fmincon','Display','none');
%     l_f_arg=fmincon(@(theta) -psi1(n,theta), arg_0, [],[], [], [], [-2 -2], [2 2],[],options);
%     %l_f_arg = fmincon(@(theta)-1* psi(n,theta), arg_0, [],[],[],options);
%     ret_val = l_f_arg;
% end

function ret_val = get_optimal_on_the_average_n(theta_vals, prob_vals,  m,sigma,c,T_star)
%  m1=5;
%  m2=5; 
A=[];
b=[];
Aeq=[c(1),c(2),c(3),c(4)];
beq=T_star;
lb=[1,1,1,1];
ub=[500,500,500,500];

des_0 = [40,30,20,10];
options = optimoptions('fmincon','Display','none');
    opt_n = fmincon(@(des)B(des,theta_vals, prob_vals, m,sigma,c,T_star), des_0,A,b,Aeq,beq,lb,ub,[],options);
    ret_val = opt_n;
end

function ret_val = B(n, theta_vals, prob_vals,  m,sigma,c,T_star)  
s = 0;
    for i = 1:length(prob_vals)
        s = s + prob_vals(i)* psi1( n, theta_vals(i,:),  m,sigma,c,T_star);
    end
    ret_val=s;
end



function y = psi1(n, theta, m,sigma,c,T_star)
for j=1:1:m
V_beta(j) =  (n)/((sigma(j)^2)*(1+((n-1)*(theta(j))-theta(j+m))));
end

y= 2/sum(V_beta);

end
