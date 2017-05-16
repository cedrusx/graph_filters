function y = FIR(co,L,x,METHOD)
% FIR filter for graph signals
% y = co(1) * x + co(2) * L * x + co(3) * L * L * x + ...

if nargin < 4 || METHOD == 0    % Centralized
    y = polynomial(co,L)*x;
elseif METHOD == 1              % Distributed
    y = co(1) * x;
    z = x;
    for t = 2:length(co)
        z = L * z;
        y = y + co(t) * z;
    end
end

% NOTICE: distributed version generally gives exactly the same output with
% centralized version