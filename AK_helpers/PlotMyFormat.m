function sp = PlotMyFormat(titleStr,xlab,ylab,zlab,xtick,ytick,ztick,varargin)
% XLabel, YLabel, ZLabel, Xticks, Yticks, Zticks
% Creates bold labels and optionally removes ticklabels

xlabS = 12;
ylabS = 12;
zlabS = 12;
titleS = 14;
fontS = 12;
if any(strcmpi(varargin,'XLabelSize'))
    xlabS = varargin{find(strcmpi(varargin,'XLabelSize'))+1};
end
if any(strcmpi(varargin,'YLabelSize'))
    ylabS = varargin{find(strcmpi(varargin,'YLabelSize'))+1};
end
if any(strcmpi(varargin,'ZLabelSize'))
    zlabS = varargin{find(strcmpi(varargin,'ZLabelSize'))+1};
end
if any(strcmpi(varargin,'TitleSize'))
    titleS = varargin{find(strcmpi(varargin,'TitleSize'))+1};
end
if any(strcmpi(varargin,'FontSize'))
    fontS = varargin{find(strcmpi(varargin,'FontSize'))+1};
end


sp = gca;

title(titleStr);
xlabel(xlab,'fontweight','bold');
ylabel(ylab,'fontweight','bold');
zlabel(zlab,'fontweight','bold');

if strcmp(xtick,'none')
    xticklabels({});
    xticks('auto');
end
if strcmp(ytick,'none')
    yticklabels({});
    yticks('auto');
end
if strcmp(ztick,'none')
    zticklabels({});
    zticks('auto');
end


%sp.XLabel.Font = 'bold';
sp.FontSize = fontS;
sp.XLabel.FontSize = xlabS;
sp.YLabel.FontSize = ylabS;
sp.ZLabel.FontSize = zlabS;
sp.Title.FontSize = titleS;




end

