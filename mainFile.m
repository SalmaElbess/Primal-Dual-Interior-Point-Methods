clc; clear;
%% Problem 1
%%from week 8 lectures
c = [-30;-20;0;0]; A = [2,1,1,0;1,3,0,1]; b = [8;8];
[x,y,s] = init_point(A,c,b);
%mehrotra
eta = 0.99;  tolerance = 0.001;
[xs_vector_mehrotra,ss_vector_mehrotra,f_vector_mehrotra] = mehrotra(x,s,y,A,c,b,eta,tolerance);
%central path fixed step
segma = 0.5; alpha = 0.9;
[xs_vector_cp_fixed,ss_vector_cp_fixed,f_vector_cp_fixed] = central_path_fixed(x,s,y,A,c,b,segma,alpha,tolerance); 
%central path adaptive step
adaptive_alpha = true; adaptive_segma = true; segma_update_parameter = 0.5;
[xs_vector_cp_adp,ss_vector_cp_adp,f_vector_cp_adp] = central_path_adaptive(x,s,y,A,c,b,segma,alpha,tolerance,adaptive_alpha,adaptive_segma,segma_update_parameter);
%matlab_built-in function
options = optimoptions('linprog','Algorithm','interior-point','Display','iter');
[x,fval,exitflag,output] = linprog(c,[],[],A,b,zeros(length(x),1),inf,options); 
%% Problem 2
% from week 8 lectures
A = [1 1 1]; c = [-1.1; -1; 0]; b = 6;
[x,y,s] = init_point(A,c,b);
%mehrotra
eta = 0.99;  tolerance = 0.001;
[xs_vector_mehrotra,ss_vector_mehrotra,f_vector_mehrotra] = mehrotra(x,s,y,A,c,b,eta,tolerance); 
%central path fixed step
segma = 0.5; alpha = 0.9;
[xs_vector_cp_fixed,ss_vector_cp_fixed,f_vector_cp_fixed] = central_path_fixed(x,s,y,A,c,b,segma,alpha,tolerance); 
%central path adaptive step
adaptive_alpha = true; adaptive_segma = true; segma_update_parameter = 0.5;
[xs_vector_cp_adp,ss_vector_cp_adp,f_vector_cp_adp] = central_path_adaptive(x,s,y,A,c,b,segma,alpha,tolerance,adaptive_alpha,adaptive_segma,segma_update_parameter); 
%matlab_built-in function
options = optimoptions('linprog','Algorithm','interior-point','Display','iter');
[x,fval,exitflag,output] = linprog(c,[],[],A,b,zeros(length(x),1),inf,options); 
%% Problem 3
%%from week 4 
c = [4;1;0;0]; A  = [3,1,0,0;4,3,-1,0;1,2,0,1]; b = [3;6;4];
[x,y,s] = init_point(A,c,b);
%mehrotra
eta = 0.99;  tolerance = 0.001;
[xs_vector_mehrotra,ss_vector_mehrotra,f_vector_mehrotra] = mehrotra(x,s,y,A,c,b,eta,tolerance); 
%central path fixed step
segma = 0.5; alpha = 0.9;
[xs_vector_cp_fixed,ss_vector_cp_fixed,f_vector_cp_fixed] = central_path_fixed(x,s,y,A,c,b,segma,alpha,tolerance); 
%central path adaptive step
adaptive_alpha = true; adaptive_segma = true; segma_update_parameter = 0.5;
[xs_vector_cp_adp,ss_vector_cp_adp,f_vector_cp_adp] = central_path_adaptive(x,s,y,A,c,b,segma,alpha,tolerance,adaptive_alpha,adaptive_segma,segma_update_parameter); 
%matlab_built-in function
options = optimoptions('linprog','Algorithm','interior-point','Display','iter');
[x,fval,exitflag,output] = linprog(c,[],[],A,b,zeros(length(x),1),inf,options); 