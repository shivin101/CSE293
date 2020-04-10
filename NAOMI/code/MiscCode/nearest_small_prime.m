function N = nearest_small_prime(N,P)

% N = nearest_small_prime(N,P)
% 
% Find the closest, larger number to N with no factors greater than P.
% 
% 2018 - Adam Charles (from code by Chris Turnes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check inputs

if isempty(P)
    P = 7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find next largest value with a good factorization

if numel(N) > 1
    for kk = 1:numel(N)
        N(kk) = nearest_small_prime(N(kk),P);
    end
else
    while max(factor(N)) > P
        N = N + 1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%