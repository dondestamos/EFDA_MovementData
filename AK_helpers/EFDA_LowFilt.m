function y = EFDA_LowFilt(x,FS,FilterParams,opts)



MaxNanFraction = opts.MaxNanInterp;
filttype = FilterParams.Type;

% simple filtering usually requires finite input, i.e. nans are not acceptable.
% if nans are at the edge of a time-series (and make a small fraction of the whole
% duration), the time-series could be truncated and still filtered, with a later replacement of the same nans. 
% But that is not reasonable for EFDA, as fdasrvf cannot handle time-series with nans. So in case
% of nans at the edges, whole time-series shall be discarded.
flag_edgenansdiscard = 1;



if strcmpi(filttype,'SavGol')
    %lpFiltF =
elseif strcmpi(filttype,'ButterLow')
    Fcutoff = FilterParams.FCutoff;
    ButterOrder = FilterParams.Order;
    [b, a] = butter(ButterOrder,Fcutoff / (FS/2));
end




%% NaN handling. Remove series which exceed a threshold of nan amount, or interpolate. Keep track of those.
Xin = x;
if ~isempty(MaxNanFraction) && nnz(isnan(x)) > 0
    inanend = zeros(1,size(x,2));
    inanbeg = ones(1,size(x,2)) * length(x) + 1;
    % Split into columns, for easier parsing
    for icol = 1:size(Xin,2)
        XcolCell{icol} = Xin(:,icol);
    end
    PostInterpNans = cell(size(Xin,2));
    for icol = 1:size(Xin,2)
        Xcol = XcolCell{icol};
        if nnz(isnan(Xcol)) > 0
            % Check if excessive - then all nans
            nanspos = find(isnan(Xcol));
            if nnz(isnan(Xcol)) > length(Xcol) * MaxNanFraction
                Xcol = nan(length(Xcol),1);
                XcolCell{icol} = Xcol;
                continue
            end
            % Check if at the edges - then truncate or discard
            if any(nanspos == length(x))
                if flag_edgenansdiscard
                    Xcol = nan(length(Xcol),1);
                    XcolCell{icol} = Xcol;
                    continue
                end
                inanbeg(1,icol) = find((isnan(Xcol(2:end,1)) & ~isnan(Xcol(1:end-1,1))),1,'last');
                Xcol(inanbeg(1,icol):end) = [];
            end
            if any(nanspos == 1)
                if flag_edgenansdiscard
                    Xcol = nan(length(Xcol),1);
                    XcolCell{icol} = Xcol;
                    continue
                end
                inanend(1,icol) = find((isnan(Xcol(1:end-1,1)) & ~isnan(Xcol(2:end,1))),1,'first');
                Xcol(1:inanend(1,icol)) = [];
            end

            % Check if any nans left - then interpolate
            nansposremain = find(isnan(Xcol));
            if isempty(nansposremain)
                XcolCell{icol} = Xcol;
                continue
            end
            % then interpolate
            irangetrunc = 1:((inanbeg(1,icol)-1) - (inanend(1,icol)+1));
            irangetruncNNans = find(~isnan(Xcol(irangetrunc)));
            irangetruncNans = find(isnan(Xcol(irangetrunc)));
            Xcol(irangetrunc) = interp1(irangetrunc(irangetruncNNans),Xcol(irangetruncNNans),irangetrunc,'spline','extrap');
            XcolCell{icol} = Xcol;
            PostInterpNans{icol} = irangetrunc(irangetruncNans);
        end
    end
end



%% Do the filtering. Replace original nans
if isempty(MaxNanFraction) || nnz(isnan(x)) == 0


    if strcmpi(filttype,'SavGol')
        %y = filtfilt(lpFiltF,x);
    elseif strcmpi(filttype,'ButterLow')
        y = filtfilt(b,a,x);
    end



else
    % Filter each column
    for icol = 1:size(Xin,2)
        if all(isfinite(XcolCell{icol}))
            if strcmpi(filttype,'SavGol')
                XColCellFilt{icol} = filtfilt(lpFiltF,XcolCell{icol});
            elseif strcmpi(filttype,'ButterLow')
                XColCellFilt{icol} = filtfilt(b,a,XcolCell{icol});
            end
        else
            XColCellFilt{icol} = nan(size(Xin,1),1);
            continue
        end
        
        % Was this column interpolated?
        if ~isempty(PostInterpNans{icol})
            temp = XColCellFilt{icol};
            temp(PostInterpNans{icol}) = nan;
            XColCellFilt{icol} = temp;
        end

        % Was this column truncated?
        if inanend(icol) + 1 > 1
            temp = XColCellFilt{icol};
            temp1 = nan(inanend(icol),1);
            temp1(end+1:end+length(temp)) = temp;
            XColCellFilt{icol} = temp1;
        end
        if inanbeg(icol) - 1 < length(Xin)
            temp = XColCellFilt{icol};
            temp(end+1:length(Xin)) = nan;%(length(Xin)-inanbeg(icol),1);
            XColCellFilt{icol} = temp;
        end

    end
    % Build columns back together
    for icol = 1:size(Xin,2)
        y(:,icol) = XColCellFilt{icol};
    end

end


end