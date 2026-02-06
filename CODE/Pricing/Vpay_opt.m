function V = Vpay_opt(y_star, zeta_e, alpha, R_fix, ~, D, D0, H, H0)
% Vpay_opt Computes optimal payer swaption price using Hagan paper approach
%   V = Vpay_opt(y_star, zeta_e, alpha, R_fix, ~, D, D0, H, H0)
%
%   INPUTS:
%     y_star  - Optimizer root 
%     zeta_e  - Integrated variance term (>0)
%     alpha   - Vector of accrual factors [alpha1,...,alphaN]
%     R_fix   - Strike swap rate
%     ~       - Unused input placeholder for market swap rates
%     D       - Vector of discount factors at payment dates
%     D0      - Initial discount factor (t=0)
%     H       - Vector of integral of mean reversion from 0 to each T
%     H0      - Integral at time 0 (usually zero)
%
%   OUTPUT:
%     V       - Payer swaption price at time 0


% Validate variance input
if zeta_e <= 0
   error('Vpay_opt: zeta_e deve essere > 0');
end

% Standard deviation of integrated process
sqz = sqrt(zeta_e);

% First term
term0 = D0 * normcdf( -y_star / sqz );

% Vector of H differences for each payment
dtH = H - H0;  % vector [H1-H0, …, Hn-H0]

arg   = (-y_star - dtH * zeta_e) / sqz;      % 1×n
coeff = alpha .* (R_fix ) .* D;           % 1×n
% Sum weighted option terms across all payments
term1 = sum( coeff .* normcdf(arg) );        % sum i=1→n
% Final term at last payment date
arg_n = (-y_star - dtH(end)*zeta_e) / sqz;
term2 = D(end) * normcdf(arg_n);

% Combine terms
V = term0 - term1 - term2;
end
