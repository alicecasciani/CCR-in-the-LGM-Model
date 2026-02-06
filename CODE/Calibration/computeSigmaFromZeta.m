function sigma = computeSigmaFromZeta(zeta, today, expiry)
% COMPUTESIGMAFROMZETA   Recover piecewise σ_i from calibrated ζ_e values
%   sigma = computeSigmaFromZeta(zeta, today, expiry)
%
%   INPUTS:
%     zeta    – 1×N vector of calibrated ζ_e at each swaption expiry (cumulative variance)
%     today   – Valuation date (datetime scalar)
%     expiry  – 1×N vector of swaption expiry dates (datetime)
%
%   OUTPUT:
%     sigma   – 1×N vector of instantaneous volatilities on each interval
%
%   The relationship is:
%     ζ_e(i) ≈ ∫_{0}^{T_i} σ(u)^2 du,
%   and so σ_i on [T_{i−1}, T_i] ≈ sqrt( (ζ_i − ζ_{i−1}) / ΔT_i ).
%   Here ΔT_i is the year‐fraction between expiries i−1 and i.

% Day-count convention
ACT_365 = 3; 
% Prepend zero to ζ so that z(1)=0
z        = [0, zeta];
% Compute year‐fractions from today to each expiry
yf1      = yearfrac(today, expiry(1), ACT_365);
% Compute year‐fractions between successive expiries
yf_rest  = yearfrac(expiry(1:end-1), expiry(2:end), ACT_365);
% Full ΔT vector
yf       = [yf1, yf_rest];
% Approximate σ on each piece: σ_i = sqrt((ζ_i − ζ_{i−1}) / ΔT_i)
z_deriv  = (z(2:end) - z(1:end-1)) ./ yf;
sigma    = sqrt(z_deriv);
end