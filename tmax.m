function t = tmax(thres, f, s, step)
% thres: threshold of normalized residual sum of squares
% f: discrete transfer function values
% s: s = 2 for IDIIR; s = 1 for FastIDIIR

if ~exist('step', 'var')
    step = 2/(max(f.^s) + min(f.^s));
end

% tset = (log(thres) - log((1-f).^2)) ./ log((1-step*f.^s).^2);
t1 = -((log(thres) - log((1-f).^2)));
t2 = -1./log((1-step*f.^s).^2);
tset = t1 .* t2;
t = max(tset);
% t = max(t2);