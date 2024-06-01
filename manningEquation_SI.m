function vf = ManningEquation_SI(n, af, pw, ss)
% This function calculates the velocity (V) in an open channel using Manning's equation.

% Inputs:
%   n: Manning's roughness coefficient (-)
%   af: flow area (m^2)
%   pw: wetted perimeter (m)
%   ss: stream slope (hf/L) (dimensionless)
%   See also: https://en.wikipedia.org/wiki/Manning_formula
% Output:
%   vf: Flow velocity (m/s)

% Check for valid input types
if ~isscalar(n) || ~isscalar(af) || ~isscalar(pw) || ~isscalar(ss)
  error('All inputs must be scalar values.');
end

if n <= 0 || af <= 0 || pw <= 0 || ss <= 0
  error('All inputs must be positive values.');
end

% Manning's equation
kappa = 1.0; % kappa = 1.0 for SI unit.
vf = kappa / n * (af / pw) ^ (2./3.) * ss ^ 0.5;
end
