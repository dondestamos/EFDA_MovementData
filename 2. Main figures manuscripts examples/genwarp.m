function [psii, gam, shvec] = genwarp(A,f,time,varargin)

% A is amplitude between 0 and 1. 0 is identity. 1 is there's at least one point
% where psi = 0 (which is technically not good).
% f is frequency of shooting vector (and of psi), see below
% Inverse flag adds pi to inverse the gamma (mirror about identity). No inverse
% means the first part (closest to zero) is convex.
% Shift is phase shift of shooting vector (and of psi), see below

% So far I've identified these, possibly useful base shapes of gamma.

% f = .5, shift = pi/2 --- no inflexion
% f = 1, shift = 0 --- no inflexion but total smaller warp (edges align to Identity)
% f = 2, shift = 0 --- all above/below identity, .5 fixed
% f = 3, shift = 0 --- all above/below identity, .33 .66 fixed
% f = 4, shift = 0 --- all above/below identity, .25 .5 .75 fixed
%...

% f = 1, shift = pi/2 --- symmetric around identity, .5 fixed
% f = 1.5 shift = pi/2 --- symmetric around identity, .33 .67 fixed
% f = 2, shift = pi/2 --- symmetric around identity, .25, .5, .75 fixed
% f = 2.5, shift = pi/2 --- symmetric around identity, .2 .4 .6 .8 fixed
% f = 3, shift = pi/2 --- symmetric around identity, 1/6 2/6 3/6 4/6 5/6 fixed
% ...


% see examples in EFDA_warps_visual_continuous_asymm('NoClose','NoSave')
% and EFDA_warps_visual_continuous_symm('NoClose','NoSave')

% Aleksei Krotov
% Northeastern University
% 2024.



if any(strcmpi(varargin,'InverseFlag'))
    InverseShift = pi * varargin{find(strcmpi(varargin,'InverseFlag'))+1};
else
    InverseShift = (randi(2) - 1) * pi; % If nothing specified, it's pseudorandom
end

if any(strcmpi(varargin,'Shift'))
    PhaseShift = varargin{find(strcmpi(varargin,'Shift'))+1};
else
    PhaseShift = 0;
end

if any(strcmpi(varargin,'PsiNoiseSigma'))
    PsiNoiseSigma = varargin{find(strcmpi(varargin,'PsiNoiseSigma'))+1};
else
    PsiNoiseSigma = 0;
end

if any(strcmpi(varargin,'PsiNoiseSigmaMulti'))
    PsiNoiseSigmaMulti = varargin{find(strcmpi(varargin,'PsiNoiseSigmaMulti'))+1};
else
    PsiNoiseSigmaMulti = 0;
end

shvec = A * sqrt(2) * atan(1/sqrt(2)) * sin(f * 2*pi .* time + InverseShift + PhaseShift);
psii = exp_map(ones(1,length(time)), shvec);

if PsiNoiseSigmaMulti ~= 0
    psii = psii .* (1 + randn(length(psii),1) .* PsiNoiseSigmaMulti);
end

if PsiNoiseSigma ~= 0
    psii = psii + randn(length(psii),1) * PsiNoiseSigma;
end

gam = cumtrapz(time,psii.^2);
gam = (gam - gam(1)) ./ (gam(end) - gam(1));


end

