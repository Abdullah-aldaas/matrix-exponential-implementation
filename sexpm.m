% Input: square matrix A
% Output: approximation of e^A using subdiagonal Padé approximant.
function r = sexpm(A)
    sigma = eigs(A,1,'lr');            % Estimate rightmost ev.
    sigma = max(sigma,0);
    A_sigma = A-sigma*eye(length(A));  % A_sigma = A-sigma*I
    no = normest(A_sigma,0.3);         % Norm estimation
    [s,k,m] =choose(no);               % Parameters dependent on the norm
    r = rkm(k,m,A_sigma/(2^s));        % Evaluation of the Padé approximant
    r = r^(2^s)*exp(sigma);            % Squaring e^A = e^sigma*e^A_sigma
end



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




% This function computes the enumerator of the (k,m)-type Padé approximant 
% Input: (k,m): Padé parameters
% Output: Polynomial p in a symbolic expression
function p = pkm(k,m)
    syms z
    p=0;
    for j=0:k
        term = factorial(k+m-j)*factorial(k)/(factorial(k+m)*factorial(k-j) ...
            *factorial(j))*z^j;
        p = p + term;
   end
end

% This function computes the denominator of the (k,m)-type Padé approximant.
% Input: (k,m): Padé parameters
% Output: Polynomial q in a symbolic expression
function q = qkm(k,m)
    syms z
    q = 0;
    for j = 0:m
        term = factorial(k+m-j)*factorial(m)/(factorial(k+m)*factorial(m-j) ...
            *factorial(j))*(-z)^j;
        q = q+ term;
    end
end





% This function computes the poles and the residues for evaluating the Padé approximant.
% Input: (k,m)
% Output: poles and residues 
function [poles,residues] = parameters(k,m)
    syms z
    p = pkm(k,m);
    q = qkm(k,m);
    c = coeffs(q,z,'All');
    poles = roots(c);
    q_prime = diff(q);
    residues = subs(p/q_prime,poles);
end



% This function computes the Padé approximant in partial fraction form.
% Input: (k,m,A)
% Output: r_km(A) in partial fraction form
function r = rkm(k,m,A)
    n = length(A);
    [poles, residues] = parameters(k,m);
    poles = double(poles);
    residues = double(residues);
    r = zeros(n,n);
    for i =1:m
        bi = poles(i);
        ai = residues(i);
        summand = ai*inv(A-bi*eye(n));
        r = r+summand;
    end
    % here is the case of superdiagonal Padé approximant, where 
    % norm of A is O(1).
    if m < k
        rrr = pkm(k,m)/(qkm(k,m));
        alpha_0 = subs(rrr,0);
        for i=1:m
            bi = poles(i);
            ai = residues(i);
            alpha_0 = alpha_0 + ai/bi;
        end
        alpha_1 = subs(rrr,1) - alpha_0;
        for i=1:m
            bi = poles(i);
            ai = residues(i);
            alpha_1 = alpha_1 - ai/(1-bi);
        end
        alpha_0 = double(alpha_0);
        alpha_1 = double(alpha_1);
        r = r + alpha_0*eye(n) + alpha_1*A;
    end
    
end
