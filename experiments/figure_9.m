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






rng(0);
n = 50;
D = -logspace(0,2,n);
x_axis = logspace(0,6,50);
s_er=[];
e_er = [];
A_norm = [];
AX_cond = [];
kexp = [];
sss = [];
counter =0;
for i = 1:50
    
    c = x_axis(i);                  % k(X)
    X = gallery('randsvd',n,c);
    X_inv = inv(X);
    A = X *diag(D) * X_inv;
    EXP = X * diag(exp(D)) * X_inv;
    [rel,est] = funm_condest1(A,@expm);
    kexp(end+1) = eps*rel;
    exp_mat = sexpm(A);
    expm_new = expm(A);
    e_norm = norm(EXP,2);
    s_er(end+1) = norm(exp_mat-EXP,2)/e_norm;
    e_er(end+1) = norm(expm_new-EXP,2)/e_norm;
    norm_of_A = norm(A,2);
    [s,k,m] = choose(norm_of_A);
    sss(end+1) = s+max(k,m);
    A_norm(end+1) = eps* norm_of_A;
    AX_cond(end+1) = eps*c*norm_of_A;
    counter = i

end





seed =[];
eee = [];

for i =1:50
    nor = 5*10^(i/5);
    [s,k,m] = choose(nor);
    s2 = find_s(nor);
    eee(end+1) = s2+13;
    seed(end+1) =i;
end









figure;
subplot(1,2,1)
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5 5 12 8]);  % 6 cm wide, 8 cm high

loglog(x_axis,e_er,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
hold on
loglog(x_axis,s_er,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
loglog(x_axis,A_norm,'b--',LineWidth=1.5)
loglog(x_axis,AX_cond,'m--',LineWidth=1.5)
loglog(x_axis,kexp,'black', LineWidth=1.5)

legend('expm','sexpm','$u\|\mathbf{A}\|$','$u\|\mathbf{A}\| \kappa_2(X)$',...
    'Interpreter', 'latex','Location','northwest','Fontsize',13)

yticks([1e-15 1e-10 1e-5 1e0])
yticklabels({'10^{-15}', '10^{-10}', '10^{-5}', '10^{0}'})
ylim([1e-15,1e0])
xlabel('$\kappa_2(X)$','Interpreter', 'latex')
ylabel('error')
pbaspect([1.5,1,1])






subplot(1,2,2)
semilogx(x_axis,sss,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
semilogx(x_axis,eee,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
legend('sexpm','expm','location','northwest', 'Fontsize', 13)
xlabel('$\kappa_2(X)$','Interpreter', 'latex')
ylabel('flops/$n^3$','Interpreter', 'latex')
pbaspect([1.5,1,1])
