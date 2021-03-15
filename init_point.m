function [x,y,s] = init_point(A,c,b)
%based on algorithm in Numerical optimization, Nocedal & Wright p.410 
x = A'*((A*A')^-1)*b;
y = ((A*A')^-1)*A*c;
s = c- A'*y;
% to ensure non-negativity
delta_x = max([0,-3*min(x)/2]);
delta_s = max([0,-3*min(s)/2]);

x = x+delta_x;
s = s+delta_s;
%avoid dissimilarity and being close to zero
delta_x_hat = x'*s/(2*sum(s));
delta_s_hat = x'*s/(2*sum(x));

x = x+delta_x_hat;
s = s+delta_s_hat;
end