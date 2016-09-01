function [ x ] = cfun( t )
% Function of precision matrix entries
% This example is a piecewise linear function
if 0.2 <= t && t <= 0.3
    x = 0.6 - 12*abs(t-0.25);
elseif 0.5 <= t && t <= 0.6
    x = -0.6 + 12*abs(t-0.55); 
elseif 0.8 <= t && t <= 0.9
    x = 0.6 - 12*abs(t-0.85);
else x = 0;
end

end

