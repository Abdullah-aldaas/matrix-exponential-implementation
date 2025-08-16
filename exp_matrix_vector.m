% Input: norm of a matrix
% Output: parameters s,(k,m)
% we assume norm is larger 1
function [s,k,m] = choose(norm)
    N = [1,200,10^4,10^6,10^9,10^11,10^12,10^14,inf];
    S = [4,4,4,3,2,2,2,1,1];
    K = [5,4,3,3,3,2,1,1,1];
    M = [4,5,4,4,4,3,2,2,2];
    for j = 2:length(N)
        if (norm >=N(j-1)) && (norm < N(j))
            s = S(j-1);
            k = K(j-1);
            m = M(j-1);
        end
    end
end



%polynomials p,q in symbolic expression
function p = pkm(k,m)
    syms z
    p=0;
    for j=0:k
        term = factorial(k+m-j)*factorial(k)/(factorial(k+m)*factorial(k-j) ...
            *factorial(j))*z^j;
        p = p + term;
   end
end

function q = qkm(k,m)
    syms z
    q = 0;
    for j = 0:m
        term = factorial(k+m-j)*factorial(m)/(factorial(k+m)*factorial(m-j) ...
            *factorial(j))*(-z)^j;
        q = q+ term;
    end
end




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



% Input: (k,m,s,A,v)
% Output: r_km(A/2^s)^{2^s}v
function rr = rkm(k,m,s,A,v)
    n = length(A);
    [poles,residues] = parameters(k,m);
    V = {v};
    for j=2:2^s+1
        r = zeros(n,1);        
        for i =1:m
            bi = poles(i);
            ai = residues(i);
            summand = (A-bi*eye(n)) \ V{j-1};
            r = r+(ai*summand);
        end

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
            r = r+alpha_0*V{j-1};
            r = r+alpha_1*A*V{j-1};
        end
        V{end+1} = r;
    end    
    rr = V{end};
end




% Input: square matrix A
% Output: approximation of e^A using subdiagonal Padé approximant.
function r = sexpmv(A,v)
    sigma = eigs(A,1,'lr');             % Estimate rightmost ev.
    sigma = max(sigma,0);
    A_sigma = A-sigma*eye(length(A));   % A_sigma = A-sigma*I
    no = normest(A_sigma,0.3);          % Norm estimation
    [s,k,m] =choose(no);                % Parameters dependent on the norm
    r = rkm(k,m,s,A_sigma/(2^s),v);     % Evaluation, 2^s linear systems
    r = r*exp(sigma);                   % e^Av ≈ e^sigma*e^A_sigma*v
end