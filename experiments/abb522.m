% This function computes the scaling paramter and the Padé parameters based on the norm of the matrix
% Input norm: norm of a matrix
% Output: parameters: 
% s: scaling parameter
% k, m: Padé parameters
% we assume norm is larger than 1
function [s,k,m] = choose(normA)
    if normA < 200
        s = 4; k = 5; m = 4;
    elseif normA < 10^4
        s = 4; k = 4; m = 5;
    elseif normA < 10^6
        s = 4; k = 3; m = 4;
    elseif normA < 10^9
        s = 3; k = 3; m = 4;
    elseif normA < 10^(11)
        s = 2; k = 3; m = 4;
    elseif normA < 10^(12)
        s = 2; k = 2; m = 3;
    elseif normA < 10^(14)
        s = 2; k = 1; m = 2;
    else 
        s = 1; k = 1; m = 2;
    end
end




function s = find_s(num)
    s = log2(num);
    s = ceil(s)+1;
end






n = 100;
rng(0);
x_axis = [];    % seed values =1:50
e_er = [];      % relative error expm
s_er = [];      % relative error sexpm
kx = [];        % u*k(X)*norm_of_A
cond_exp = [];  % u*k_exp
sss = [];
eee = [];
v = rand(n,1);
for seed=1:50
    A = generate_random(seed,n);
    [X,D] = eig(A);
    X_inv = inv(X);
    EXP = X * diag(exp(diag(D))) * X_inv *v;
    [rel,est] = funm_condest1(A,@expm);
    cond_exp(end+1) = eps*rel;
    norm_of_A = norm(A,2);
    kx(end+1) = eps*norm_of_A*cond(X);
    exp_mat = sexpmv(A,v);
    expm_new = expmv(A,v);
    e_norm = norm(EXP,2)*norm(v,2);
    s_er(end+1) = norm(exp_mat-EXP,2)/e_norm;
    e_er(end+1) = norm(expm_new-EXP,2)/e_norm;
    [s,k,m] = choose(norm_of_A);
    s2 = find_s(norm_of_A);
    sss(end+1) = 2^s*m;
    eee(end+1) = 2^(s2)*13;
    x_axis(end+1) = seed;
    seed
end






figure;
subplot(1,2,1)
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5 5 12 8]);  % 6 cm wide, 8 cm high

semilogy(x_axis,e_er,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
hold on
semilogy(x_axis,s_er,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
semilogy(x_axis,kx,'m--',LineWidth=1.5)
semilogy(x_axis,cond_exp,'black', LineWidth=1.5)

legend('expmv','sexpmv','$u\|\mathbf{A}\| \kappa_2(X)$',...
    'Interpreter', 'latex','Location','northwest','Fontsize',13)

yticks([1e-15 1e-10 1e-5 1e0])
yticklabels({'10^{-15}', '10^{-10}', '10^{-5}', '10^{0}'})
ylim([1e-15,1e0])
xlabel('seed')
ylabel('error')
pbaspect([1.5,1,1])




subplot(1,2,2)
semilogy(x_axis,sss,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
semilogy(x_axis,eee,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
legend('sexpmv','expmv','location','northwest', 'Fontsize', 13)
xlabel('seed')
ylabel('linear systems')
pbaspect([1.5,1,1])