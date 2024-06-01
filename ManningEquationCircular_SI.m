function vf = ManningEquationCircular_SI(n, diameter, depth, S)
% This function calculates the velocity (V) in an open channel using Manning's equation.

% Inputs:
%   n: Manning's roughness coefficient 
%   A: flow area (m^2)
%   diameter: diameter of circular tube
%   depth: depth of water, which should be between 0 and diameter.
%          If depth > diameter, this function assumes it is diameter, 
%          and prints a warning message. 
%   S: stream slope (hf/L) (dimensionless)
%   See also: https://en.wikipedia.org/wiki/Manning_formula
% Output:
%   V: Flow velocity (m/s)

% Check for valid input types
if ~isscalar(n) || ~isscalar(diameter) || ~isscalar(depth) || ~isscalar(S)
  error('All inputs must be scalar values.');
end

if n <= 0 || diameter <= 0 || depth <= 0 || S <= 0
  error('All inputs must be positive values.');
end

% Manning's equation
if depth > diameter
    fprintf('# Warning from ManningEquationCircular_SI.m: diameter=%f, depth=%f, but depth should be <= diameter.\n', diameter, depth);
end
%  calculate rk (hudraulic radius), af (flow area), and pw (wetted
%  perimeter)
[rh, af, pw] = circularFlowHydraulicRadius(0.5 * diameter, depth); 
kappa = 1.0; % kappa = 1.0 for SI unit.
vf = kappa / n * (af / pw) ^ (2./3.) * S ^ 0.5;
end
