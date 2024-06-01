function ht_t_predicted = ...
    predictHeight_v6_qin_qout(ht_init, rain_mm_sub, diameter, ...
    hb, kappa, gamma, ckns, cqin)
% This function calculates change of water level (height) (per minute).
% Given the future rain data, the initial water level, and optimized 
% parameters, this function predicts future water level. 
%
% Inputs:
%   ht_init: The initial value of height (water level), which is the
%            It is a scalar, not an array.
%   rain_mm_sub: selected part of rain history (mm per .. hour?) 
%            (time series by minute) 
%            The rain_mm_sub here does not need to be the entire rain_mm. 
%            It only needs to be a section of rain data that you want to
%            predict the height (water level).
%            The dimension of rain_mm_sub is (nt, 1). 
%            (nt - 1) is the length of water level you want to predict.
%   diameter: diameter of main tube (horizontal, slightly slope)
%   hb, kappa, gamma, ckns, cqin:
%            Optimized parameters.
% Output: 
%   ht_t_predicted: the predicted height (water level). 
%                   ht_t_predicted(1) is not predicted. It is ht_init.
%                   Size of the returned ht_t_predicted is the same size 
%                   with rain_mm. 
    % number of steps to predict (precisely speaking, the index 1 is not
    % predicted. It is the initial value.) 
    nt = size(rain_mm_sub(:), 1);
    % run explicit calculation (numerical integral)
    % numerical integral
    ht_explicit = zeros(nt,1);
	ht_explicit(1) = ht_init; 
    for i = 2: nt
		dhdt_i = dhdt_v5_qin_qout(...
			ht_explicit(i - 1), ...
			.5 * (rain_mm_sub(i - 1) + rain_mm_sub(i)), diameter, ...
            hb, kappa, gamma, ckns, cqin);	
		ht_explicit(i) = ht_explicit(i - 1) ...
    		+ dhdt_i * 1.0;
    end 
  	ht_t_predicted = ht_explicit(:);
end
