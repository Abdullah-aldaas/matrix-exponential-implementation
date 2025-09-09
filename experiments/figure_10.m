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



n = 50;
rng(0);
X = gallery('randsvd',n,10);
X_inv = inv(X);
spread=[];
s_er=[];
e_er = [];
A_norm = [];
kexp = [];
v= rand(n,1);
sss = [];
eee = [];
for i = 1:50   
    spread(end+1) = 10^(i/5);
    D = -logspace(0,i/5,n);
    D = D + randn/20*j*D;
    A = X * diag(D) * X_inv;
    EXP = X*diag(exp(D))*X_inv *v;
    [c,est] = funm_condest1(A,@expm);
    kexp(end+1) = eps*c;
    exp_mat = sexpmv(A,v);
    expm_new = expmv(A,v);
    e_norm = norm(EXP,2)*norm(v,2);
    s_er(end+1) = norm(exp_mat-EXP,2) / e_norm;
    e_er(end+1) = norm(expm_new-EXP,2) / e_norm;
    norm_of_A = norm(A,2);
    A_norm(end+1) = eps* norm_of_A;
    [s,k,m] = choose(norm_of_A);
    s2 = find_s(norm_of_A);
    sss(end+1) = 2^s*m;
    eee(end+1) = 2^(s2)*13;
    counter = i

end



figure;
subplot(1,2,1)
set(gcf, 'Position', [5 5 12 8]);  % 6 cm wide, 8 cm high

loglog(spread,s_er,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
loglog(spread,e_er,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
loglog(spread,A_norm,'b--',LineWidth=1.5)
loglog(spread,kexp,'black', LineWidth=1.5)
legend('sexpmv','expmv','$u\|\mathbf{A}\|$', 'Interpreter', 'latex','Location','southeast','Fontsize',13)

yticks([1e-15 1e-10 1e-5 1e0])
yticklabels({'10^{-15}', '10^{-10}', '10^{-5}', '10^{0}'})
ylim([1e-15,1e0])
xlabel('spread')
ylabel('error')
pbaspect([1.5,1,1])




subplot(1,2,2)
loglog(spread,sss,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
loglog(spread,eee,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
legend('sexpmv','expmv','location','northwest', 'Fontsize', 13)
xlabel('spread')
ylabel('linear systems')
pbaspect([1.5,1,1])
