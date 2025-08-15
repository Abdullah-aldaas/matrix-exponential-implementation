% This function generates matrix A with dimension n, that is :
% essentially nonnegative.
% diagonal dominant.
% all eigenvalues of A are <=0. 
% ein eigenvalue is close to 0.
% norm of A is ca. 10^(seed/5).
function A = generate_mat(seed,n)
    rng(seed);
    A = rand(n,n);
    A = 10^(seed/5)*A/norm(A);
    A = abs(A);
    
    % this is to ensure there is an eigenvalue close to 0.
    for j =2:n
        A(1,j) = abs(10^(-9)*rand());
    end
    
    % make diagonal elements small enough.
    for i =1:n
        s = 0;
        for j = 1:n
            if i ~= j
                s= s+A(i,j) ;
            end
        end
        A(i,i) = -s-rand();
    end
end


