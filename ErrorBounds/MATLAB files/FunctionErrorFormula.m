
% n: number of derivatives of a function used
% Mvalues: maximum magnitude of (n+1) derivatives of function in the range of
% interest; for example, if n = 1, then we are interested in the maximum
% value of the n+1 = 2nd derivative of the function

errors_n = nan(Nstates,1);

nn = n:-1:n-1;
coeff = ts.^(nn+1)./ factorial(nn+1);


% From data
Mvalues = max_local_errors(1:Nstates)'./coeff;
% Mvalues = max(max_local_errors(1:6)*2/0.01^2, max_local_errors(14:26)*1/0.01);
% Analytical
% Mvalues = [16.52; 15.86; 16.06; 12.11; 12.12; 12.11; 12.11; 42.37; 60.05; 68.74; 269.42; 60.40; 43.64];

% Mvalues = [
for states = 1 : Nstates
        errors_n(states,1) = Mvalues(states) * (m*ts)^(nn(states)+1)/factorial(nn(states)+1);
end
errors_n

