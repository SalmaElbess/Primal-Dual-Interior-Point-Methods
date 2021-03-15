function [xs_vector,ss_vector,f_vector] = central_path_adaptive(x,s,y,A,c,b,segma,alpha,tolerance,adaptive_alpha,adaptive_segma,segma_update_parameter) 
%calculations
[m,n] = size(A);
%----------------
n_iterations = 0;
xs_vector = x;
ss_vector = s;
f_vector = c'*x;
while sum(x.*s) >= tolerance
    %updating parameters
    u = sum(x.*s)/n;
    T = segma*u;
    n_iterations = n_iterations+1;
    %diag X S
    X = diag(x);
    S = diag(s);
    %jaccobian
    j = [zeros(n) transpose(A) eye(n);
        A        zeros(m,m)   zeros(m,n);
        S        zeros(n,m)   X];
    %right hand side
    rc = transpose(A)*y + s - c;
    rb = A*x - b;
    %solving
    delta = j\[-1*rc; -1*rb ;-1*(x.*s)+T];
    x_delta = delta(1:n);
    y_delta = delta(n+1:n+m);
    s_delta = delta(n+m+1:end);
    %updating alpha using mehrotra hurestics   
    if adaptive_alpha
    alpha_pri = min([1; min(-1*x(x_delta<0)./x_delta(x_delta<0))]);
    alpha_dual = min([1; min(-1*s(s_delta<0)./s_delta(s_delta<0))]);
    %updating points
    x = x + alpha_pri*delta(1:n);
    y = y + alpha_dual*delta(n+1:n+m);
    s = s + alpha_dual*delta(n+m+1:end);
    else
    x = x + alpha*delta(1:n);
    y = y + alpha*delta(n+1:n+m);
    s = s + alpha*delta(n+m+1:end);
    end
    
    %updating segma 
    if adaptive_segma
        segma = segma*segma_update_parameter;
    end
    
    %alpha 
    xs_vector = [xs_vector , x];
    ss_vector = [ss_vector,s];
    f_vector = [f_vector , c'*x];
end
figure
plot(xs_vector(1,:),xs_vector(2,:),'-*')
xlabel('x1'); ylabel('x2')
title('Central Path with adaptive step size and centering parameter')

figure
plot(xs_vector(1,:).*ss_vector(1,:),xs_vector(2,:).*ss_vector(2,:),'-*')
xlabel('x1*s1'); ylabel('x2*s2')
title('Central Path with adaptive step size and centering parameter')

figure
plot(0:1:n_iterations,f_vector,'-*')
xlabel('iteration'); ylabel('objective function')
title('Central Path with adaptive step size and centering parameter')

end