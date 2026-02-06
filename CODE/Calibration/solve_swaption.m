function sol = solve_swaption(alpha, R_fix, S, D, D0, H, H0, mkt_price)
% SOLVE_SWAPTION    Solve for (y*, ζ_e) via nonlinear least squares
%   sol = solve_swaption(alpha, R_fix, S, D, D0, H, H0, mkt_price)
%
%   We want to find v = [ y_star; ζ_e ] such that two equations are satisfied:
%     F1(v) = 0  (present‐value‐matching)
%     F2(v) = 0  (swaption‐price‐matching)
%
%   INPUTS:
%     alpha     – 1×n vector of accrual fractions for each fixed‐leg payment
%     R_fix     – Scalar: fixed swap rate for this swaption
%     S         – (Not used here; could be placeholder for forward rate info)
%     D         – 1×n vector of P(0, t_k) for payment dates t1…tn
%     D0        – Scalar: P(0, expiry)
%     H         – 1×n vector of H(τ_k) for payment dates
%     H0        – Scalar: H(τ_0) at expiry
%     mkt_price – Scalar: observed market swaption price
%
%   OUTPUT:
%     sol       – 2×1 vector [y_star; ζ_e] that minimizes ||F(v)||^2


% Define the residual‐vector function F = [F1; F2]
fun_resid = @(v) equations(v, alpha, R_fix, S, D, D0, H, H0, mkt_price);

% Initial guess: y_star ≈ 0.01, ζ_e ≈ 0.01
v0 = [0.01; 0.01];

% Lower/upper bounds: y_star free, ζ_e ≥ 0
lb = [-Inf; 0];
ub = [ Inf; Inf];

% Opzioni per lsqnonlin
 options = optimoptions('lsqnonlin', ...
    'Display', 'off', ...                      % no output
    'Algorithm', 'trust-region-reflective', ...% supporta bounds
    'FunctionTolerance', 1e-10, ...
    'StepTolerance',     1e-10, ...
    'MaxIterations',    100);


% Call lsqnonlin to find optimal v
[v_opt, resnorm, res, exitflag, output] = ...
    lsqnonlin(fun_resid, v0, lb, ub, options);

if exitflag <= 0
    warning('lsqnonlin non ha convergito (exitflag = %d)', exitflag);
end

sol = v_opt;
end