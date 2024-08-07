function q = f_to_srvf(f,time,spl)
% F_TO_SRVF Convert function to Square-Root Velocity Function
% -------------------------------------------------------------------------
% Convert to SRSF
%
% Usage: q = f_to_srvf(f,time)
%
% This function converts functions to srsf
%
% Input:
% f: matrix of functions
% time: vector of time samples
% spl: use b-spline computation (default: true)
%
% Output:
% q: matrix of SRSFs
if nargin < 3
    spl = true;
end

binsize = mean(diff(time));
[M, N] = size(f);

if spl
    fy = zeros(M,N);
    for ii = 1:N
        
        
        
        % AKR: What if I use mydiff instead with central differences and interp at the end? - MUCH BETTER
        % Nans are not expected to make it to EFDA. 
%         if nnz(isnan(f(:,ii))) > 0.1 * M
%             y = Bspline(f(:,ii),3);
%             ydiff = diff(y);
%             fy(:,ii) = ydiff(1:length(time))/binsize;
%             %f(:,ii) = y(1:length(time)); % This is not further used, why is it here??? //Alex
%         else
%             fy(:,ii) = mydiff(time,f(:,ii));
%         end
        fy(:,ii) = mydiff_finite(time,f(:,ii));
        

    end
    q = fy./sqrt(abs(fy)+eps);
else % Don't even use this one.
    if size(f,2)>1
        [~,fy] = gradient(f,binsize,binsize);
    else
        fy = gradient(f,binsize);
    end
    q = fy./sqrt(abs(fy)+eps);
end

end




function y = mydiff_finite(t,x)
% Differentiation via central differences

assert(iscolumn(t));
assert(size(x,1) == size(t,1));
ndof = size(x,2);

ind0 = isnan(x(:,1));   

dt = (t(3:end) - t(1:end-2)) ./ 2;
dx = (x(3:end,:) - x(1:end-2,:)) ./ 2;
yy = dx ./ dt;
yy(2:end+1,:) = yy;
yy(1,:) = (x(2,:) - x(1,:)) ./ (t(2) - t(1));
yy(end+1,:) = (x(end,:) - x(end-1,:)) ./ (t(end) - t(end-1));
t1 = linspace(t(1),t(end),length(t))';

y = nan(size(x,1),size(x,2));
for idof = 1:ndof
    y(:,idof) = interp1(t,yy(:,idof),t1,'makima',nan); % ,'extrap' % Extrap sometimes shoots the result in infinity.
end

y(ind0,:) = nan;

end

