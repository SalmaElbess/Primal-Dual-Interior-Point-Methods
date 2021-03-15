function [xs_vector,ss_vector,f_vector] = central_path_fixed(x,s,y,A,c,b,segma,alpha,tolerance) 
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
    %updating points
    x = x + alpha*delta(1:n);
    y = y + alpha*delta(n+1:n+m);
    s = s + alpha*delta(n+m+1:end);
    %alpha 
    xs_vector = [xs_vector , x];
    ss_vector = [ss_vector,s];
    f_vector = [f_vector , c'*x];
end
figure
plot(xs_vector(1,:),xs_vector(2,:),'-*')
xlabel('x1'); ylabel('x2')
title('Central Path with fixed step size and centering parameter ')
figure
plot(xs_vector(1,:).*ss_vector(1,:),xs_vector(2,:).*ss_vector(2,:),'-*')
xlabel('x1*s1'); ylabel('x2*s2')
title('Central Path with fixed step size and centering parameter')
figure
plot(0:1:n_iterations,f_vector,'-*')
xlabel('iteration'); ylabel('objective function')
title('Central Path with fixed step size and centering parameter')

end