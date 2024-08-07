function [sp, os] = TiledAxesMy(curFig,NRows,NCols,optionsIn,varargin)
% os is short for options
if nargin < 4
    optionsIn = [];
end


if nargin < 1
    NRows = 1;
    NCols = 1;
end

% Create default properties
os.NRows = NRows;
os.NCols = NCols;
os.FigPosNorm = [0.01 0.01 0.99 0.99];
os.SpacingsOut = [40 50 50 50]; % [Right Top Left Bottom] in pixels. Default ~60
os.SpacingsOutStr = 'Right Top Left Bottom in px'; %. Default ~60
os.SpacingsBetween = [50 50];
os.SpacingsBetweenStr = 'Vert Horiz';
os.YLabels = 'First'; % May change to "All", "First" or to "FirstEvery2nd" or "BotLeftCorn"
os.XLabels = 'First'; % May change to "All", "First" or to "FirstEvery2nd" or "BotLeftCorn"
os.XTickMarks = 'All'; % Same
os.YTickMarks = 'All'; % Same
os.NPlotsHoriz = NRows * NCols; %how many from (Nrows*Ncols) should be added? The last row(s) can be incomplete.
% os.NPlotsVert = NRows * NCols; %how many from (Nrows*Ncols) should be added? The last column(s) can be incomplete.
% The last option is not yet supported!

if nargin < 1
    sp = os;
    return
end
optFlds = fieldnames(os);
if ~isempty(optionsIn) && isstruct(optionsIn)
    optInFlds = fieldnames(optionsIn);
    for ifld = 1:length(optInFlds)
        if ismember(optInFlds{ifld},optFlds)
            os.(optInFlds{ifld}) = optionsIn.(optInFlds{ifld});
        else
            disp(sprintf('<strong>TiledAxesMy: Unknown option field supplied: %s</strong>',optInFlds{ifld}));
        end
    end
end

if isempty(curFig)
    curFig = figure;
end
figure(curFig);

% Compute Dimensions
SR = os.SpacingsOut(1); % right
ST = os.SpacingsOut(2); % top
SL = os.SpacingsOut(3); % left
SB = os.SpacingsOut(4); % bottom
SV = os.SpacingsBetween(1); %vert
SH = os.SpacingsBetween(2); %horiz
DimsOuterPx = floor(os.FigPosNorm(3:4) .* curFig.Position(3:4));
AxW = floor((DimsOuterPx(1) - SR - SL - (NCols-1) .* SH) ./ NCols) ./ curFig.Position(3); % modify if YLabels not all.
AxH = floor((DimsOuterPx(2) - ST - SB - (NRows-1) .* SV) ./ NRows) ./ curFig.Position(4); % modify if YLabels not all.
SpaceOut([1 3]) = os.SpacingsOut([1 3]) ./ curFig.Position(3);
SpaceOut([2 4]) = os.SpacingsOut([2 4]) ./ curFig.Position(4);
SpaceBet = os.SpacingsBetween ./ curFig.Position([4 3]);

% Add more data to options, if needs to be used for other plots in the figure...
os.ComputedDimensionsUnits = 'Figure-normalized';
os.AxesWidth = AxW;
os.AxesHeight = AxH;




% Where labels
xlabCol = 1:NCols;
xlabRow = 1:NRows;
if strcmpi(os.XLabels,'First')
    xlabRow = NRows;
elseif strcmpi(os.XLabels,'FirstEvery2nd')
    xlabRow = NRows;
    xlabCol = 1:2:NCols;
elseif strcmpi(os.XLabels,'BotLeftCorn')
    xlabRow = NRows;
    xlabCol = 1;
elseif strcmpi(os.XLabels,'None')
    xlabRow = [];
    xlabCol = [];
end
ylabCol = 1:NCols;
ylabRow = 1:NRows;
if strcmpi(os.YLabels,'First')
    ylabCol = 1;
elseif strcmpi(os.YLabels,'FirstEvery2nd')
    ylabCol = 1;
    ylabRow = NRows:(-2):1;
elseif strcmpi(os.YLabels,'BotLeftCorn')
    ylabCol = 1;
    ylabRow = NRows;
elseif strcmpi(os.YLabels,'None')
    ylabCol = [];
    ylabRow = [];
end

% Where ticklabels
xtickCol = 1:NCols;
xtickRow = 1:NRows;
if strcmpi(os.XTickMarks,'First')
    xtickRow = NRows;
elseif strcmpi(os.XTickMarks,'FirstEvery2nd')
    xtickRow = NRows;
    xtickCol = 1:2:NCols;
elseif strcmpi(os.XTickMarks,'BotLeftCorn')
    xtickRow = NRows;
    xtickCol = 1;
elseif strcmpi(os.XTickMarks,'None')
    xtickRow = [];
    xtickCol = [];
end
ytickCol = 1:NCols;
ytickRow = 1:NRows;
if strcmpi(os.YTickMarks,'First')
    ytickCol = 1;
elseif strcmpi(os.YTickMarks,'FirstEvery2nd')
    ytickCol = 1;
    ytickRow = NRows:(-2):1;
elseif strcmpi(os.YTickMarks,'BotLeftCorn')
    ytickCol = 1;
    ytickRow = NRows;
elseif strcmpi(os.YTickMarks,'None')
    ytickCol = [];
    ytickRow = [];
end


% Create plots
for irow = 1:NRows
    for icol = 1:NCols
        isp = (irow-1) * NCols + icol;
        PosX = SpaceOut(3) + (icol-1) * (AxW + SpaceBet(2)) + os.FigPosNorm(1);
        PosY = SpaceOut(4) + (NRows-irow) * (AxH + SpaceBet(1)) + os.FigPosNorm(2);
        sp(isp) = subplot('Position',[PosX PosY AxW AxH]);
        os.AxPositions(:,isp) = [PosX PosY];
        hold on;
        
        if ismember(icol,xlabCol) && ismember(irow,xlabRow)
            xlabel('Label','fontweight','bold');
        end
        if ismember(icol,ylabCol) && ismember(irow,ylabRow)
            ylabel('Label','fontweight','bold');
        end
        if ~(ismember(icol,xtickCol) && ismember(irow,xtickRow))
            xticklabels({});
            xticks('auto');
        end
        if ~(ismember(icol,ytickCol) && ismember(irow,ytickRow))
            yticklabels({});
            yticks('auto');
        end
        
        sp(isp).YLabel.FontWeight = 'bold';
        sp(isp).XLabel.FontWeight = 'bold';
        
        if os.NPlotsHoriz < NRows * NCols && isp == os.NPlotsHoriz
            break
        end
    end
end

end