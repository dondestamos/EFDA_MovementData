function [y,Status] = mydiff(t,x,varargin)
% Numerical differentiation first order
% Central differences with edge-padding by one value (i1 = i2, iend = iend-1)
% Returns uniform sampling even if input had non-uniform sampling

% Takes vector of time (may be not equidistantly spaced) Nx1 or 1xN,
% Takes vector of values or a matrix of values N x M
% Preferably, N is a column for both. Automatic reshaping works if M < N.

% ! Add 2nd and 3rd derivatives?

% Aleksei Krotov
% Northeastern University, 2023

Status.Val = 0;
Status.Msg = '';
warning('off','MATLAB:interp1:NaNstrip');
s = size(x);

flag_filter = 0;
if any(strcmpi(varargin,'LPFilter'))
	flag_filter = 1;
    FS = varargin{find(strcmpi(varargin,'LPFilter'))+1};
    Fcut = varargin{find(strcmpi(varargin,'LPFilter'))+2};
end


if any(strcmpi(varargin,'Method'))
	InterpMethod = varargin{find(strcmpi(varargin,'Method'))+1};
else
	InterpMethod = 'makima';  %Could be changed (e.g. to makima) to avoid undulations
end

if any(strcmpi(varargin,'ExtrapArg'))
	ExtrapArg = varargin{find(strcmpi(varargin,'ExtrapArg'))+1};
else
	ExtrapArg = nan; % By defrault, no extrapolation.
end

flag_separateColumns = 0;
if any(strcmpi(varargin,'SeparateColumnsForNanPurpose'))
	flag_separateColumns = 1;
end

flag_Speed = 0;
if any(strcmpi(varargin,'Speed'))
	flag_Speed = 1;
end

if any(strcmpi(varargin,'DerivMethod'))
	DerivMethod = varargin{find(strcmpi(varargin,'DerivMethod'))+1};
else
	DerivMethod = 'numeric_spline';  % Could be changed (to 'analytic_spline') to ensure continuity. 
end


if any(strcmpi(varargin,'NoReshape'))
    ndof = s(2);
    %fprintf('<strong>Make sure mydiff arguments are time-columns!</strong>');
else
    
%Reshapes automatically.
ndof = min(s);
% columns expected
if s(1) < s(2) %I WANT IT TO BE COLUMN(S)
    x = x';
end
y = nan(size(x,1),size(x,2));
if flag_Speed
    y = nan(size(x,1),1);
end

if (nnz(isnan(x)) > 0.8 * numel (x))
    Status.Val = 1;
    Status.Msg = 'mydiff: More than 80%% nans';
    return
end
st = size(t);
if st(1) < st(2)
    t = t';
end
assert(max(st) == max(s),'Time and data lengths are not equal!');
end


if flag_separateColumns && ndof > 1
    ind0 = zeros(s);
    ind0 = ind0 == 0;
    for idof = 1:ndof
        xcurr = x(:,idof);
        % Remember nan locations
        ind0(:,idof) = isnan(xcurr(:,1));   
    end
else
    ind0 = isnan(x(:,1));   
end


% Remember nan locations
ind0 = isnan(x(:,1));   
ind = isnan(diff(x,1,1));


% Also try upsampling x2, taking 1+i/2 : end-i/2, then re-interpolating and resampling. 
% Check results on a sine signal of a specific frequency. Try different sampling
% frequencies, to see how these methods correspond to Nyquist. 
% Compare them with forward and backward differences.

%%%%%% Central differences then stretch - produce higher error, esp. when sampling frequency is
%%%%%% low compared with the harmonics frequencies
% dt = (t(3:end) - t(1:end-2)) ./ 2;
% dx = (x(3:end) - x(1:end-2)) ./ 2;
% yy = dx ./ dt;
% t1 = linspace(t(1),t(end),length(t)-2)';

%%%%%% Forward differences, then stretch - produce half the error of central diff.
% dt = t(2:end) - t(1:end-1); 
% dx = x(2:end,:) - x(1:end-1,:);
% yy = dx ./ dt;
% t1 = linspace(t(1),t(end),length(t)-1)';

%%%%%% Central differences without stretch, edges extrapolated from second-order diff,
%%%%%% produce about 1/10 of the error of central diff with stretch.
% dt = (t(3:end) - t(1:end-2)) ./ 2;
% dx = (x(3:end,:) - x(1:end-2,:)) ./ 2;
% yy = dx ./ dt;
% ddx = x(3:end,:) + x(1:end-2,:) - 2*x(2:end-1,:);
% ddt = (t(3:end) - t(1:end-2)) ./ 2;
% yyy = ddx ./ ddt.^2;
% yy(2:end+1,:) = yy;
% yy(1,:) = yy(2,:) - dt(1) .* yyy(1,:);
% yy(end+1,:) = yy(end,:) + dt(end) .* yyy(end,:);
% t1 = linspace(t(1),t(end),length(t))';

%%%%%% Central differences without stretch, Edges extrapolated from forw and backw differences
%%%%%% produce about 1/20 of the error of central diff with stretch
dt = (t(3:end) - t(1:end-2)) ./ 2;
dx = (x(3:end,:) - x(1:end-2,:)) ./ 2;
yy = dx ./ dt;
yy(2:end+1,:) = yy;
yy(1,:) = (x(2,:) - x(1,:)) ./ (t(2) - t(1));
yy(end+1,:) = (x(end,:) - x(end-1,:)) ./ (t(end) - t(end-1));
t1 = linspace(t(1),t(end),length(t))';


if strcmpi(DerivMethod,'numeric_spline')
for idof = 1:ndof
    if all(isnan(yy(:,idof)))
        continue
    end
    y(:,idof) = interp1(t,yy(:,idof),t1,InterpMethod,ExtrapArg); % ,'extrap' % Extrap sometimes shoots the result in infinity.
end
elseif strcmpi(DerivMethod,'analytic_spline') % To be continued maybe. Still not understand completely.
    error('Analytic spline not supported');
    for idof = 1:ndof
        ppspl = interp1My(t,x(:,idof),InterpMethod,'pp'); % ,'extrap' % I think this one is not correct
        %y(:,idof) = interp1(t1,yy(:,idof),t,InterpMethod); % ,'extrap' % I think this one is not correct
    end
end


if flag_filter
    y = LowFilt(y,Fcut - Fcut/10,Fcut + Fcut/10,'FS',FS,'FilterType','Butter4');
end

warning('on','MATLAB:interp1:NaNstrip');

if flag_separateColumns && ndof > 1
    for idof = 1:ndof
        y(ind0(:,idof),idof) = nan;
    end
else
    y(ind0,:) = nan;
end

if flag_Speed
    y = sqrt(sum(y.^2,2));
end


