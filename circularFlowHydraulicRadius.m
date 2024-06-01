function [rh, af, pw] = circularFlowHydraulicRadius(r, depth)
% This function calculates the flow area of a circular pipe (A) and 
% the wetted perimeter (p)

% Inputs:
%   r: radius of pipe
%   depth: depth of water (0:no water. r:half full. 2*r: full)
%   rh: hydraulic radius = (flow area) / (wetted perimeter)
%   af: flow area
%   pw: wetted perimeter

    % Check for valid input types
    if ~isscalar(r) || ~isscalar(depth)
      error('All inputs must be scalar values.');
    end
    if r < 0
      error('Radius must be a positive value.');
    end
    if depth < 0
      error('Depth must be a positive value (between 0 and 2r).');
    end
    if depth > (2 * r)
      error('Depth cannot be greater than 2r (between 0 and 2r).');
    end

    % theta
    if depth < r
      h = max(0.00000001, depth);
    else
      h = 2 * r - depth;
    end
    theta = 2 * acos((r - h) / r);
    % circular segment area
    K = 0.5 * r * r * (theta - sin(theta));
    % flow area af, wetted perimeter pw, and hydraulic radius rh
    if depth < r
      af = K;
      pw = r * theta;
    else
      af = pi * r * r - K;
      pw = 2 * pi * r - r * theta;  
    end
    % rh
    rh = af / pw; 
end