function y = LowFilt(x, Fpass, Fstop,varargin)



x = cast(x,'double');

if any(strcmpi(varargin,'FilterType'))
    filttype = varargin{find(strcmpi(varargin,'FilterType'))+1};
else
    filttype = 'FIR';
end
if any(strcmpi(varargin,'FS'))
    Fs = varargin{find(strcmpi(varargin,'FS'))+1};
else
    Fs = 500;
end

% Allow to specify the fraction of nans (between 0 and 1) that will be temporarily
% interpolated, then replaced with nans again.
if any(strcmpi(varargin,'MaxNanNTol'))
    MaxNanFraction = varargin{find(strcmpi(varargin,'MaxNanNTol'))+1};
    assert(MaxNanFraction >= 0 && MaxNanFraction <=1);
else
    MaxNanFraction = 0.25;
end




% nancheck = isnan(x);
% if sum(sum(nancheck)) > 0
%     y = zeros(size(x,1),size(x,2));
%     return
% end


% 
% %Fpass = 20;   % Passband Frequency
% %Fstop = 21;   % Stopband Frequency
% Apass = 1;    % Passband Ripple (dB)
% Astop = 60;   % Stopband Attenuation (dB)
% Fs    = 500;  % Sampling Frequency
% 
% h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, 1, 60, 500);
% 
% Hd = design(h, 'butter', ...
%     'MatchExactly', 'stopband', ...
%     'SOSScaleNorm', 'Linf');

%%
% Fpass = 39;
% Fstop = 41;
% Fs = 500;

if strcmp(filttype,'FIR')
    lpFiltF = designfilt('lowpassfir','PassbandFrequency',Fpass,...
        'StopbandFrequency',Fstop,'PassbandRipple',0.1, ...
        'StopbandAttenuation',60,'DesignMethod','kaiserwin',...
        'SampleRate',Fs);
elseif strcmp(filttype,'IIR')
    lpFiltI = designfilt('lowpassiir','FilterOrder',4, ...
         'PassbandFrequency',(Fpass+Fstop)/2,'PassbandRipple',0.1, ...
         'SampleRate',Fs);
elseif strcmp(filttype,'Butter4')
    [b, a] = butter(4,(Fpass+Fstop)/2 / (Fs/2));     
elseif strcmp(filttype,'Butter2')
    [b, a] = butter(2,(Fpass+Fstop)/2 / (Fs/2));     
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
            % Check if at the edges - then truncate
            if any(nanspos == length(x))
                inanbeg(1,icol) = find((isnan(Xcol(2:end,1)) & ~isnan(Xcol(1:end-1,1))),1,'last');
                Xcol(inanbeg(1,icol):end) = [];
            end
            if any(nanspos == 1)
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
            
   
% fvtool(lpFiltF);

%%


% lpFiltI = designfilt('lowpassiir','PassbandFrequency',Fpass, ...
%     'StopbandFrequency',Fstop,'PassbandRipple',0.4, ...
%     'StopbandAttenuation',60,'DesignMethod','ellip',...
%     'SampleRate',Fs);

%fvtool(lpFiltI);

%%

% for i = 1:numpass
%     y = filtfilt(lpFiltF,x);


if isempty(MaxNanFraction) || nnz(isnan(x)) == 0
    
    if strcmp(filttype,'FIR')
        y = filtfilt(lpFiltF,x);
    elseif strcmp(filttype,'IIR')
        y = filtfilt(lpFiltI,x);
    elseif strcmp(filttype,'Butter4')
        y = filtfilt(b,a,x);
    elseif strcmp(filttype,'Butter2')
        y = filtfilt(b,a,x);
    end


else
    % Filter each column
    for icol = 1:size(Xin,2)
        if all(isfinite(XcolCell{icol}))
            if strcmp(filttype,'FIR')
                XColCellFilt{icol} = filtfilt(lpFiltF,XcolCell{icol});
            elseif strcmp(filttype,'IIR')
                XColCellFilt{icol} = filtfilt(lpFiltI,XcolCell{icol});
            elseif strcmp(filttype,'Butter4')
                XColCellFilt{icol} = filtfilt(b,a,XcolCell{icol});
            elseif strcmp(filttype,'Butter2')
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
    
    

