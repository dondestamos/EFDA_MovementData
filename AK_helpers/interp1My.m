function Vout = interp1My(X,V,Xq,method,extrapval)
%INTERP1 1-D interpolation (table lookup)
%
%   Vq = INTERP1(X,V,Xq) interpolates to find Vq, the values of the
%   underlying function V=F(X) at the query points Xq. 
%
%   X must be a vector. The length of X is equal to N.
%   If V is a vector, V must have length N, and Vq is the same size as Xq.
%   If V is an array of size [N,D1,D2,...,Dk], then the interpolation is
%   performed for each D1-by-D2-by-...-Dk value in V(i,:,:,...,:). If Xq
%   is a vector of length M, then Vq has size [M,D1,D2,...,Dk]. If Xq is 
%   an array of size [M1,M2,...,Mj], then Vq is of size
%   [M1,M2,...,Mj,D1,D2,...,Dk].
%
%   Vq = INTERP1(V,Xq) assumes X = 1:N, where N is LENGTH(V)
%   for vector V or SIZE(V,1) for array V.
%
%   Interpolation is the same operation as "table lookup".  Described in
%   "table lookup" terms, the "table" is [X,V] and INTERP1 "looks-up"
%   the elements of Xq in X, and, based upon their location, returns
%   values Vq interpolated within the elements of V.
%
%   Vq = INTERP1(X,V,Xq,METHOD) specifies the interpolation method.
%   The available methods are:
%
%     'linear'   - (default) linear interpolation
%     'nearest'  - nearest neighbor interpolation
%     'next'     - next neighbor interpolation
%     'previous' - previous neighbor interpolation
%     'spline'   - piecewise cubic spline interpolation (SPLINE)
%     'pchip'    - shape-preserving piecewise cubic interpolation
%     'cubic'    - cubic convolution interpolation for uniformly-spaced
%                  data. This method does not extrapolate and falls back to
%                  'spline' interpolation for irregularly-spaced data.
%            NOTE: 'cubic' changed in R2020b to perform cubic convolution.
%                  In previous releases, 'cubic' was the same as 'pchip'.
%     'v5cubic'  - same as 'cubic'
%     'makima'   - modified Akima cubic interpolation
%
%   Vq = INTERP1(X,V,Xq,METHOD,'extrap') uses the interpolation algorithm
%   specified by METHOD to perform extrapolation for elements of Xq outside
%   the interval spanned by X.
%
%   Vq = INTERP1(X,V,Xq,METHOD,EXTRAPVAL) replaces the values outside of
%   the interval spanned by X with EXTRAPVAL.  NaN and 0 are often used for
%   EXTRAPVAL.  The default extrapolation behavior with four input
%   arguments is 'extrap' for 'spline', 'pchip' and 'makima', and
%   EXTRAPVAL = NaN (NaN+NaN*1i for complex values) for the other methods.
%
%   PP = INTERP1(X,V,METHOD,'pp') is not recommended. Use
%   griddedInterpolant instead.
%   PP = INTERP1(X,V,METHOD,'pp') uses the interpolation algorithm
%   specified by METHOD to generate the ppform (piecewise polynomial form)
%   of V. The method may be any of the above METHOD except for 'v5cubic'
%   and 'makima'. PP may then be evaluated via PPVAL. PPVAL(PP,Xq) is the
%   same as INTERP1(X,V,Xq,METHOD,'extrap').
%
%   For example, generate a coarse sine curve and interpolate over a
%   finer abscissa:
%       X = 0:10; V = sin(X); Xq = 0:.25:10;
%       Vq = interp1(X,V,Xq); plot(X,V,'o',Xq,Vq,':.')
%
%   For a multi-dimensional example, we construct a table of functional
%   values:
%       X = [1:10]'; V = [ X.^2, X.^3, X.^4 ];
%       Xq = [ 1.5, 1.75; 7.5, 7.75]; Vq = interp1(X,V,Xq);
%
%   creates 2-by-2 matrices of interpolated function values, one matrix for
%   each of the 3 functions. Vq will be of size 2-by-2-by-3.
%
%   Class support for inputs X, V, Xq, EXTRAPVAL:
%      float: double, single
%
%   See also INTERPFT, SPLINE, PCHIP, INTERP2, INTERP3, INTERPN, PPVAL,
%            griddedInterpolant, scatteredInterpolant.

%   Copyright 1984-2020 The MathWorks, Inc.

% Parse interpolation method, extrapolation value, number of consecutive
% numeric data arguments, and whether we are in the 'pp' case or not.
pp = [];
ndataarg = nargin; % Number of X,V,Xq args. Init to nargin and reduce.

if nargin == 2 && isfloat(V)
    % INTERP1(V,Xq)
    method = 'linear';
    extrapval = 'none';
elseif nargin == 3 && isfloat(Xq)
    % INTERP1(X,V,Xq)
    if isScalarTextArg(V)
        error(message('MATLAB:interp1:nargin'));
    end
    method = 'linear';
    extrapval = 'none';
else
    % INTERP1(X,V,Xq,___)
    % INTERP1(V,Xq,___)
    if nargin == 5
        last = extrapval;
        penultimate = method;
    elseif nargin == 4
        last = method;
        penultimate = Xq;
    elseif nargin == 3
        last = Xq;
        penultimate = V;
    else
        last = V;
        penultimate = X;
    end
    [X,V,method,extrapval,ndataarg,pp] = ...
        parseinputs(X,V,penultimate,last,nargin,ndataarg,pp);
end



%%%% MY CHECHUPS
flag_nan = 0;
if all(isnan(V),'all')
    flag_nan = 1;
end

if all(isnan(Xq),'all')
    flag_nan = 1;
end

if all(isnan(X),'all')
    flag_nan = 1;
end

NNNans = nnz(~isnan(V(:)));
if NNNans < 3
    flag_nan = 1;
end

if flag_nan
    %dim = length(size(V));
    sizes = size(V);
    sizes(1) = 0;
    Vout = nan(sizes);
    Vout = repmat(Vout,length(Xq));
    return
end

% Make sure the inputs are unique
%[Xu,ic,iv] = unique(X);



warning('off','MATLAB:interp1:NaNstrip');
if strcmpi(method,'linear') && strcmpi(extrapval,'none')
    Vout = interp1(X,V,Xq);
elseif strcmpi(extrapval,'none')
    Vout = interp1(X,V,Xq,method);
else
    Vout = interp1(X,V,Xq,method,extrapval);
end
    warning('on','MATLAB:interp1:NaNstrip');
 

end



%%%%

%-------------------------------------------------------------------------%
%     'nearest'  - nearest neighbor interpolation
%     'next'     - next neighbor interpolation
%     'previous' - previous neighbor interpolation
%     'linear'   - linear interpolation
%     'spline'   - piecewise cubic spline interpolation (SPLINE)
%     'pchip'    - shape-preserving piecewise cubic interpolation
%     'cubic'    - same as 'pchip'
%     'v5cubic'  - the cubic interpolation from MATLAB 5
%     'makima'   - modified Akima cubic interpolation
function methodname = sanitycheckmethod(method_in)
method = char(method_in); % string support
if isempty(method)
    methodname = 'linear';
else
    if method(1) == '*'
        method(1) = [];
    end
    if matlab.internal.math.partialMatch(method, 'linear')
        methodname = 'linear';
    elseif matlab.internal.math.partialMatch(method, 'makima')
        methodname = 'makima';
    elseif matlab.internal.math.partialMatch(method, 'spline')
        methodname = 'spline';
    elseif matlab.internal.math.partialMatch(method, 'next', 3)
        methodname = 'next';
    elseif matlab.internal.math.partialMatch(method, 'nearest')
        methodname = 'nearest';
    elseif matlab.internal.math.partialMatch(method, 'pchip')
        methodname = 'pchip';
    elseif matlab.internal.math.partialMatch(method, 'previous', 2)
        methodname = 'previous';
    elseif matlab.internal.math.partialMatch(method, 'cubic')
        methodname = 'cubic';
    elseif matlab.internal.math.partialMatch(method, 'pp', 2) % legacy behavior
        methodname = 'pchip';
    elseif matlab.internal.math.partialMatch(method, 'v5cubic')
        methodname = 'cubic';     
    else
        error(message('MATLAB:interp1:InvalidMethod'));
    end
end
end

%-------------------------------------------------------------------------%
function pp = ppinterp(X,V, orig_size_v, method)
%PPINTERP ppform interpretation.
n = size(V,1);
ds = 1;
prodDs = 1;
if ~isvector(V)
    ds = orig_size_v(2:end);
    prodDs = size(V,2);
end

switch method(1)
    case 'n' % next and nearest
        if strcmpi(method(3), 'x')
            error(message('MATLAB:interp1:ppGriddedInterpolantNext'));
        else
            breaks = [X(1); (X(1:end-1)+X(2:end))/2; X(end)].';
            coefs = V.';
            pp = mkpp(breaks,coefs,ds);
        end
    case 'l' % linear
        breaks = X.';
        page1 = (diff(V)./repmat(diff(X),[1, prodDs])).';
        page2 = (reshape(V(1:end-1,:),[n-1, prodDs])).';
        coefs = cat(3,page1,page2);
        pp = mkpp(breaks,coefs,ds);
    case 'p' % previous, pchip and cubic
        if strcmpi(method(2), 'r')
			error(message('MATLAB:interp1:ppGriddedInterpolantPrevious'));
        else
            pp = pchip(X.',reshape(V.',[ds, n]));
        end
    case 's' % spline
        pp = spline(X.',reshape(V.',[ds, n]));
    case 'c' % v5cubic
        b = diff(X);
        if norm(diff(b),Inf) <= eps(norm(X,Inf))
            % data are equally spaced
            a = repmat(b,[1 prodDs]).';
            yReorg = [3*V(1,:)-3*V(2,:)+V(3,:); ...
                V; ...
                3*V(n,:)-3*V(n-1,:)+V(n-2,:)];
            y1 = yReorg(1:end-3,:).';
            y2 = yReorg(2:end-2,:).';
            y3 = yReorg(3:end-1,:).';
            y4 = yReorg(4:end,:).';
            breaks = X.';
            page1 = (-y1+3*y2-3*y3+y4)./(2*a.^3);
            page2 = (2*y1-5*y2+4*y3-y4)./(2*a.^2);
            page3 = (-y1+y3)./(2*a);
            page4 = y2;
            coefs = cat(3,page1,page2,page3,page4);
            pp = mkpp(breaks,coefs,ds);
        else
            % data are not equally spaced
            pp = spline(X.',reshape(V.',[ds, n]));
        end
end

% Even if method is 'spline' or 'pchip', we still need to record that the
% input data V was oriented according to INTERP1's rules.
% Thus PPVAL will return Vq oriented according to INTERP1's rules and
% Vq = INTERP1(X,Y,Xq,METHOD) will be the same as
% Vq = PPVAL(INTERP1(X,Y,METHOD,'pp'),Xq)
pp.orient = 'first';

end % PPINTERP
%-------------------------------------------------------------------------%
function [X,V,method,extrapval,ndataarg,pp] = ...
    parseinputs(X,V,penultimate,last,ninputs,ndataarg,pp)
if isScalarTextArg(last)
    if strcmp(last,'pp')
        if (ninputs ~= 4)
            error(message('MATLAB:interp1:ppOutput'))
        end
        method = sanitycheckmethod(penultimate);
        if strcmp(method,'makima')
            error(message('MATLAB:interp1:ppAkima'));
        end
        extrapval = 'none';
        [X,V,orig_size_v] = reshapeAndSortXandV(X,V);
        % Use griddedInterpolant constructor for remaining error checks
        griddedInterpolant(X,V(:,1));
        pp = ppinterp(X, V, orig_size_v, method);
        
    elseif strcmp(last,'extrap')
        if (ninputs ~= 4 && ninputs ~= 5)
            error(message('MATLAB:interp1:nargin'));
        end
        if ~isempty(penultimate) && ~isScalarTextArg(penultimate)
            error(message('MATLAB:interp1:ExtrapNoMethod'));
        end
        method = sanitycheckmethod(penultimate);
        ndataarg = ninputs-2;
        extrapval = 'extrap';
        if strcmp(method,'cubic')
            extrapval = 'none';
            warning(message('MATLAB:interp1:NoExtrapForV5cubic'))
        end
    else
        if isScalarTextArg(penultimate)
            error(message('MATLAB:interp1:InvalidSpecPPExtrap'))
        end
        method = sanitycheckmethod(last);
        extrapval = 'none';
        needextrapval = ~any(strcmp(method,{'spline','pchip','makima'}));
        if ~needextrapval
            extrapval = 'extrap';
        end
        ndataarg = ninputs-1;
    end
elseif isscalar(last)
    if isScalarTextArg(penultimate)
        extrapval = last;
        ndataarg = ninputs-2;
        method = sanitycheckmethod(penultimate);
    elseif isempty(penultimate) && (ninputs == 4 || ninputs == 5)
        % default method via []
        method = 'linear';
        extrapval = last;
        ndataarg = ninputs-2;
    else
        method = 'linear';
        extrapval = last;
    end
elseif isempty(last)
    % This is potentially ambiguous, the assumed intent is case I
    % I)    X,  V, []  Empty query
    % II)   V, [], []  Empty query and empty method
    % III)  V, Xq, []  Empty method
    method = 'linear';
    extrapval = 'none';
    if ninputs ~= 3
        ndataarg = ninputs-1;
    end
else
    method = 'linear';
    extrapval = 'none';
end
end
%-------------------------------------------------------------------------%
function [x,V,orig_size_v] = reshapeAndSortXandV(x,V)
% Reshape and sort x and V. Minimal error checking. The rest of the error
% checking is done later in the griddedInterpolant constructor.
if ~isfloat(x)
    error(message('MATLAB:interp1:Xnumeric'));
end
if ~isvector(x) % also catches empty x
    error(message('MATLAB:interp1:Xvector'));
end
[V,orig_size_v] = reshapeValuesV(V);
% Reshape x and V
x = x(:);
if numel(x) ~= size(V,1)
    if isvector(V)
        error(message('MATLAB:interp1:YVectorInvalidNumRows'))
    else
        error(message('MATLAB:interp1:YInvalidNumRows'));
    end
end
% We can now safely index into V
if ~issorted(x)
    [x,idx] = sort(x);
    V = V(idx,:);
end
end
%-------------------------------------------------------------------------%
function [V, orig_size_v] = reshapeValuesV(V)
% Reshapes V into a matrix so that we can interpolate down each column
if ~isfloat(V)
    error(message('MATLAB:interp1:NonFloatValues'));
end
orig_size_v = size(V);
if isvector(V)
    V = V(:);
elseif ~ismatrix(V)
    nrows = orig_size_v(1);
    ncols = prod(orig_size_v(2:end));
    V = reshape(V,[nrows ncols]);
end
end
%-------------------------------------------------------------------------%

function tf = isScalarTextArg(s)
tf = ischar(s) || (isstring(s) && isscalar(s));
end