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





n = 100;
X_axis = [];    % seed values =1:50
e_er = [];      % relative error expm
s_er = [];      % relative error sexpm
kx = [];        % u*k(X)*norm_of_A
cond_exp = [];  % u*k_exp
p = [];
na = [];
EE = [];
for seed=1:50
    A = generate_random(seed,n);
    [X,D] = eig(A);
    X_inv = inv(X);
    EXP = X * diag(exp(diag(D))) * X_inv;
    exp_mat = sexpm(A);
    expm_new = expm(A);
    [c,est] = funm_condest1(A,@expm);
    e_norm = norm(EXP,2);
    cx = cond(X);

    EE(end+1) = e_norm;
    p(end+1) = cx;
    s_er(end+1) = norm(EXP-exp_mat,2) / e_norm;
    e_er(end+1) = norm(EXP-expm_new,2) / e_norm;
    kx(end+1) = eps*cx*norm(A,2);
    cond_exp(end+1) = eps*c; 
    X_axis(end+1) = seed;
    na(end+1) = norm(A,2);
    seed
    
end






function s = find_s(num)
    s = log2(num);
    s = ceil(s)+1;
end



seed =[];
sss = [];
eee = [];

for i =1:50
    nor = 10^(i/5);
    [s,k,m] = choose(nor);
    s2 = find_s(nor);
    sss(end+1) = s+max(k,m);
    eee(end+1) = s2+13;
    seed(end+1) =i;
end







figure;
%set(gcf, 'Units', 'centimeters');
%set(gcf, 'Position', [5 5 12 8]);  % 6 cm wide, 8 cm high
subplot(1,2,1);
%ylim([1e-15 1e-5]); % Optional: restrict visible range
set(gcf, 'Position', [5 5 12 8]);  % 6 cm wide, 8 cm high




semilogy(X_axis,s_er,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
semilogy(X_axis,e_er,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
semilogy(X_axis,kx,'b--')
semilogy(X_axis, cond_exp,'black', LineWidth=1.5)
legend('sexpm','expm','$u\kappa_2(X)\|\mathbf{A}\|$', 'Interpreter', 'latex','Location','southeast', 'fontsize', 13)
yticks([1e-15 1e-10 1e-5 1e0])
yticklabels({'10^{-15}', '10^{-10}', '10^{-5}', '10^{0}'})
ylim([1e-15,1e0])
xlabel("seed")
ylabel("error")
pbaspect([1.5 1 1])



subplot(1,2,2)
plot(seed,sss,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
plot(seed,eee,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
legend('sexpm','expm','location','northwest', 'fontsize', 13)
xlabel('seed')
ylabel('flops/$n^3$','Interpreter', 'latex')
pbaspect([1.5,1,1])