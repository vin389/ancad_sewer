function ht_t_predicted = ...
    predictHeight_v3_qin_qout(rain_mm_sub, ht_init, ...
    diameter, gamma, cqin, lambda, hb, ckns)
% This function calculates change of water level (height) (per minute).
% Given the future rain data, the initial water level, and optimized 
% parameters, this function predicts future water level. 
%
% Inputs:
%   rain_mm_sub: selected part of rain history (mm per .. hour?) 
%            (time series by minute) 
%            The rain_mm_sub here does not need to be the entire rain_mm. 
%            It only needs to be a section of rain data that you want to
%            predict the height (water level).
%            The dimension of rain_mm_sub is (nt, 1). 
%            (nt - 1) is the length of water level you want to predict.
%   ht_init: The initial value of height (water level), which is the
%            It is a scalar, not an array.
%   diameter, gamma, cqin, lambda, hb, ckns:
%            Optimized parameters.
% Output: 
%   ht_t_predicted: the predicted height (water level). 
%                   ht_t_predicted(1) is not predicted. It is ht_init.
%                   Size of the returned ht_t_predicted is the same size 
%                   with rain_mm. 
    % number of steps to predict (precisely speaking, the index 1 is not
    % predicted. It is the initial value.) 
    nt = size(rain_mm_sub(:), 1);
    % tspan of prediction. The first element must refer to the initial
    % value. 
    tspan = [1 nt]; 
    % predict by using ode23
    [t_ode23, ht_ode23] = ode23(@(t, ht_t) ...
                dhdt_v3_qin_qout_ode23(...
                    t, ht_t, rain_mm_sub, ...
                    diameter, gamma, cqin, lambda, hb, ckns), ...
                    tspan, ht_init); 
    
	% if ode23 is okay, run interpolation
    if (sum(isnan(ht_ode23)) <= 0) 
        % interpolation 
        ht_t_predicted = interp1(t_ode23, ht_ode23, tspan);
    end
    % if ode23 fails, or interpolation fails, 
    % run explicit calculation (numerical integral)
    if (sum(isnan(ht_ode23)) > 0 || size(ht_t_predicted(:),1) ~= nt) 
        % numerical integral
        fprintf('# Warning: Ode23 failed (got nan results) in this case.\n');
        ht_explicit = zeros(nt,1);
		ht_explicit(1) = ht_init; 
        for i = 2: nt
			dhdt_i = dhdt_v3_qin_qout(...
				ht_explicit(i - 1), ...
				.5 * (rain_mm_sub(i - 1) + rain_mm_sub(i)), ...
				diameter, gamma, cqin, lambda, hb, ckns);	
			ht_explicit(i) = ht_explicit(i - 1) ...
				+ dhdt_i * 1.0;
        end 
    	ht_t_predicted = ht_explicit;
    end
end
