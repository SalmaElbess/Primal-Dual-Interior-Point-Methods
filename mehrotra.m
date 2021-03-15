function [xs_vector,ss_vector,f_vector] = mehrotra(x,s,y,A,c,b,eta,tolerance) 
%calculations
[m,n] = size(A);
%----------------
n_iterations = 0;
xs_vector = x;
ss_vector = s;
f_vector = c'*x;
while x'*s >= tolerance
    mu = x'*s/n;
    n_iterations = n_iterations+1;
    % affine scaling
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
    %solving affine scaling step
    delta_aff = j\[-1*rc; -1*rb ;-1*(x.*s)];
    delta_aff_x = delta_aff(1:n);
    delta_aff_y = delta_aff(n+1:n+m);
    delta_aff_s = delta_aff(n+m+1:end);
    % update segma
    alpha_pri_aff = min([1; min(-1*x(delta_aff_x<0)./delta_aff_x(delta_aff_x<0))]);
    alpha_dual_aff = min([1; min(-1*s(delta_aff_s<0)./delta_aff_s(delta_aff_s<0))]);
    mu_affine = sum((x + alpha_pri_aff*delta_aff_x).*(s + alpha_dual_aff*delta_aff_s))/n;
    segma = (mu_affine / mu)^3;
    % corrector step
    delta = j\[-1*rc; -1*rb ;-1*(x.*s) - (delta_aff_x.*delta_aff_s) + segma*mu];
    delta_x = delta(1:n);
    delta_y = delta(n+1:n+m);
    delta_s = delta(n+m+1:end);
    %get alpha corrector
    alpha_pri_k_max = min(-1*x(delta_x<0)./delta_x(delta_x<0));
    alpha_dual_k_max = min(-1*s(delta_s<0)./delta_s(delta_s<0));
    alpha_pri_k = min([1; eta*alpha_pri_k_max]);
    alpha_dual_k = min([1; eta*alpha_dual_k_max]);
    x = x + alpha_pri_k*delta_x;
    y = y + alpha_dual_k*delta_y;
    s = s + alpha_dual_k*delta_s;
    xs_vector = [xs_vector , x];
    ss_vector = [ss_vector,s];
    f_vector = [f_vector , c'*x];
end
figure
plot(xs_vector(1,:),xs_vector(2,:),'-*')
xlabel('x1'); ylabel('x2')
title('Mehrotra')
figure
plot(xs_vector(1,:).*ss_vector(1,:),xs_vector(2,:).*ss_vector(2,:),'-*')
xlabel('x1*s1'); ylabel('x2*s2')
title('Mehrotra')
figure
plot(0:1:n_iterations,f_vector,'-*')
xlabel('iteration'); ylabel('objective function')
title('Mehrotra')

end