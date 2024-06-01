function [rhs, afs, pws] = circularFlowHydraulicRadius_vec(r, depths)
% This function calculates the flow area of a circular pipe (A) and 
% the wetted perimeter (p)

% Inputs:
%   r: radius of pipe
%   depths: a vector of depths of water (0:no water. r:half full. 2*r: full)
%   rh: hydraulic radius = (flow area) / (wetted perimeter)
%   af: flow area
%   pw: wetted perimeter
    % theta
    h = max(0.00000001, depths);
    h = min(2. * r, h);
    % Case 1: depth < r.  Case 2: depth >= r.
    use1 = (depths < r);
    use2 = ~use1;
    % calculate h for both cases (h1 and h2)
    h1 = h;         % for h < r, small h
    h2 = 2 * r - h; % for h >= r, big h
    % calculate theta and K for both cases
    theta1 = 2 * acos((r - h1) / r);
    theta2 = 2 * acos((r - h2) / r);
    K1 = 0.5 * r * r * (theta1 - sin(theta1));
    K2 = 0.5 * r * r * (theta2 - sin(theta2));
    af1 = K1;
    pw1 = r * theta1;
    af2 = pi * r * r - K2;
    pw2 = 2 * pi * r - r * theta2;
    % af, pw, rh
    afs = use1 .* af1 + use2 .* af2;
    pws = use1 .* pw1 + use2 .* pw2;
    rhs = afs ./ pws; 
end


