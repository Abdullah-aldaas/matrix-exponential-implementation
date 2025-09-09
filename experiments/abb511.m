
% This function computes the parameters s,k,m based on the norm of A:
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
sss = [];
counter =0;
for i = 1:50    
    
    spread(end+1) = 10^(i/5);
    D = -logspace(0,i/5,n);
    %D(end+1) = 0;
    D = D + randn/20*j*D;
    A = X_inv * diag(D) * X;
    nor = norm(A,2);
    [s,k,m] = choose(nor);
    sss(end+1) = s+max(k,m);
    %A_s = vpa(A);
    %EXP = expm(A_s);
    EXP = X_inv*diag(exp(D))*X;
    [c,est] = funm_condest1(A,@expm);
    kexp(end+1) = eps*c;
    exp_mat = sexpm(A);
    expm_new = expm(A);
    e_norm = norm(EXP,2);
    s_er(end+1) = norm(exp_mat-EXP,2) / e_norm;
    e_er(end+1) = norm(expm_new-EXP,2) / e_norm;
    A_norm(end+1) = eps* norm(A,2);
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
%set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5 5 12 8]);  % 6 cm wide, 8 cm high

loglog(spread,e_er,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
hold on
loglog(spread,s_er,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
loglog(spread,A_norm,'b--',LineWidth=1.5)
loglog(spread,kexp,'black', LineWidth=1.5)
legend('expm','sexpm','$u\|\mathbf{A}\|$', 'Interpreter', 'latex','Location','southeast','fontsize',13)

yticks([1e-15 1e-10 1e-5 1e0])
yticklabels({'10^{-15}', '10^{-10}', '10^{-5}', '10^{0}'})
ylim([1e-15,1e0])
xlabel('spread')
ylabel('error')
pbaspect([1.5 1 1])








subplot(1,2,2)
semilogx(spread,sss,'o',MarkerSize=10, LineWidth=2, MarkerEdgeColor='r')
hold on
semilogx(spread,eee,'s', Markersize= 10, linewidth=2, ...
    MarkerEdgeColor='g')
legend('sexpm','expm','location','northwest', 'fontsize', 13)
xlabel('spread')
ylabel('flops/$n^3$','Interpreter', 'latex')
pbaspect([1.5,1,1])

