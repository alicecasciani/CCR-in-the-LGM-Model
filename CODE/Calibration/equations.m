function F = equations(v, alpha, R_fix, S, D, D0, H, H0, mkt_price)
% EQUATIONS   Defines the two equations (F1, F2) for solving (y_star, ζ_e)
%   F = equations(v, alpha, R_fix, S, D, D0, H, H0, mkt_price)
%
%   INPUTS:
%     v         – 2×1 vector [y_star; ζ_e]
%     alpha     – 1×n accrual factors for each fixed leg payment
%     R_fix     – Scalar: fixed swap rate for the swap
%     S         – (unused here)
%     D         – 1×n vector of discount factors P(0, t_k) for t1…tn
%     D0        – Scalar: P(0, expiry)
%     H         – 1×n vector H(τ_k) for t1…tn
%     H0        – Scalar: H(τ_0) at expiry
%     mkt_price – Scalar: observed market swaption price
%
%   OUTPUT:
%     F         – 2×1 vector [F1; F2] to be driven to zero

y_star = v(1);
zeta_e = v(2);
F = zeros(2,1);

n = length(alpha);

% --- Equation 1
term1 = 0;
for i = 1:n
    term1 = term1 + alpha(i) * (R_fix ) * D(i) * exp( -(H(i)-H0).*y_star - 0.5*(H(i)-H0)^2*zeta_e);
end
term1 = term1 + D(n) * exp(-(H(n)-H0).*y_star - 0.5*(H(n)-H0)^2.*zeta_e );

F(1,1) = term1 - D0;

% --- Equation 2:
term2 = 0;
for i = 1:n
    arg = (-y_star - (H(i)-H0)*zeta_e)/sqrt(zeta_e);
    term2 = term2 + alpha(i)*(R_fix ) * D(i) * normcdf(arg);
end
arg_n = (-y_star - (H(n)-H0)*zeta_e)/sqrt(zeta_e);
term2 = term2 + D(n) * normcdf(arg_n);

arg_0 = -y_star / sqrt(zeta_e);
term2 = D0 * normcdf(arg_0) - term2-mkt_price;

F(2,1) = term2;

end