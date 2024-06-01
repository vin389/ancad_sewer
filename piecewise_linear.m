% It is like a Heaviside step function, except that
%   1. Rather than jumping at x=0, user can determine the jumping location
%   of x. 
%   2. Instead of a sudden jump, it linearly increases from x=x1 (y=0) to
%   x=x2 (y=1). 
function px = piecewise_linear(x1, x2, x)
	px = x(:);
	for i = 1: size(px, 1)
		if px(i) < x1
			px(i) = 0.;
		elseif(px(i)> x2)
			px(i) = 1.;
		else
			px(i) = (px(i) - x1) / (x2 - x1);
		end
	end
end

