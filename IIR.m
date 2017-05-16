function y = IIR(nu,de,L,x,METHOD,tmax,FORM,step,frq)
% IIR filter for graph signals
% y = inv(de(1)*I + de(2)*L + de(3)*L*L + ...)
%      * (nu(1)*I + nu(2)*L + de(3)*L*L + ...) * x
% METHOD
%   0 = (default) Centralized IIR 
%   1 = Iterative Distributed IIR (IDIIR)
%   2 = Fast Iterative Distributed IIR (FastIDIIR)
%   3 = Distributed IIR by Chebyshev approximation (CDIIR)
% tmax
%  for IDIIR/FastIDIIR
%   if tmax >= 1, it's used as number of iterations in each IDIIR/FastIDIIR
%   if tmax < 1, it's used as the target of NRSS to calculate tmax
%  for CDIIR
%   maximal order of Chebyshev polynomials
% FORM (only for IDIIR/FastIDIIR)
%   0 = (default) Direct form
%   1 = Fully cascade form
%   2 = Fully parallel form
% step (only for IDIIR/FastIDIIR) 
%   if step > 0, it's used as the step length
%   if step = 0 or unassigned, the optimal step length is used
% frq
%   freqencies to determine the optimal step length in IDIIR/FastIDIIR
%   freqencies to determine the range of Chebyshev approximation in CDIIR

DOMAIN = 'R';
% DOMAIN = 'C': Decomposition of cascade/parallel form will be over C
% DOMAIN = 'R': Decomposition of cascade/parallel form will be over R

if nargin < 5 || METHOD == 0    % Centralized IIR
    y = polynomial(nu,L)*(polynomial(de,L)\x);
    return
elseif METHOD == 3  % CDIIR
    if ~exist('frq','var')
        frq = eig(L);
    end
    y = ChebFilter(@(x)rational(nu,de,x),L,x,tmax,[min(frq),max(frq)]);
%     % Step 1. Calculation of coefficients
%     c = chebcoeffs(chebfun(@(x)rational(nu,de,x),[min(frq),max(frq)],tmax));
%     % Step 2. Distributed filtering
%     y = chebpolyvalm(c,L,[min(frq),max(frq)]) * x;
    return
end

% IDIIR/FastIDIIR
if nargin < 7
    FORM = 0;
end
if nargin < 8
    step = 0;
end
if nargin < 9
    frq = eig(L);
end
if FORM == 0
    y = FIR(nu,L,x);
    y = DIIR(de,L,y,METHOD,tmax,step,frq);
elseif FORM == 1
    y = FIR(nu,L,x);
    poles = sort(roots(de(end:-1:1)));
    if DOMAIN == 'C'
        for k = 1:length(poles)
            y = DIIR([1 -1/poles(k)],L,y,METHOD,tmax,step,frq);
        end
    else
        k = 1;
        while k <= length(poles)
            if k < length(poles) && poles(k) == poles(k+1)'
                sde = [1, (-poles(k)-poles(k+1))/(poles(k)*poles(k+1)), 1/(poles(k)*poles(k+1))];
                k = k + 2;
            else
                sde = [1, -1/poles(k)];
                k = k + 1;
            end
            y = DIIR(sde,L,y,METHOD,tmax,step,frq);
        end
    end
    y = real(y / de(1));
elseif FORM == 2
    [snu,sde,tail] = partialfraction(nu,de,DOMAIN);
    y = zeros(size(x));
    for k = 1:length(snu)
        z = DIIR(sde{k}/sde{k}(1),L,x,METHOD,tmax,step,frq);
        y = y + FIR(snu{k}/sde{k}(1),L,z);
    end
    if ~isempty(tail)
        y = y + FIR(tail,L,x);
    end
    y = real(y);
end

function y = DIIR(de,L,x,METHOD,tm,step,frq)
B = polynomial(de,L);
if METHOD == 0  % Centralized IIR
    y = B\x;
    return;
end
if nargin < 6 || step == 0
    if nargin < 7
        hf = polynomial(de,eig(L));
    else
        hf = polynomial(de,frq);
    end
    step = 2 / (min(hf.^(3-METHOD)) + max(hf.^(3-METHOD)));
end
if tm < 1
    tm = tmax(tm,polynomial(de,eig(L)),3-METHOD,step);
end
y = x;
% ========================= Naive Approach ================================
% if METHOD == 1  % IDIIR
%     for t = 1:tm
%         y = y - step * B * (B * y - x);
%     end
% elseif METHOD == 2 % FastIDIIR
% %     if min(polynomial(de,eig(L))) < -1e-4
% %         warning(['Denominator of frequency response is not positive definite. ' ...
% %             'FastIDIIR may not converge! Minimum: %g'], min(h));
% %     end
%     for t = 1:tm
%         y = y - step * (B * y - x);
%     end
% end
% ========================= Optimized Approach ============================
if METHOD == 1  % IDIIR
    P = eye(length(x)) - step * B * B;
    b = step * B * x;
else    % FastIDIIR
    P = eye(length(x)) - step * B;
    b = step * x;
end
for t = 1:tm
    y = P * y + b;
end
% =========================================================================