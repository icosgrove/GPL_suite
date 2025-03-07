function setLabel(figPARAMS,label,n,fig,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets the x or y labels in requested locations
% Written: Ian 10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes1 = gca;
switch label
    case 'x'
        switch figPARAMS.labelLoc
            case 'All'
                xlabel(figPARAMS.xlabel)
                axes1.FontSize = figPARAMS.xFontSize;
                if figPARAMS.xlabelBold == 1
                    axes1.FontWeight = 'bold';
                end
                set(axes1,'FontWeight','normal')
            case 'Outside'
                if (figPARAMS.col(n)==1) || (figPARAMS.row(n)==figPARAMS.subplotSize(1)) % proceed
                    xlabel(figPARAMS.xlabel)
                    axes1.FontSize = figPARAMS.xFontSize;
                    if figPARAMS.xlabelBold == 1
                        axes1.FontWeight = 'bold';
                    end
                end
            case 'One'
                if (figPARAMS.row(n)==1) && (figPARAMS.col(n)==1)
                    han = axes(fig,'visible','off'); 
                    han.XLabel.Visible='on';
                    xlabel(han,figPARAMS.xlabel);
                    axes1.FontSize = figPARAMS.xFontSize;
                    if figPARAMS.xlabelBold == 1
                        axes1.FontWeight = 'bold';
                    end
                    set(axes1,'FontWeight','normal')
                end
        end
    case 'y1'
            switch figPARAMS.labelLoc
                case 'All'
                    ylabel(figPARAMS.y1label)
                    axes1.FontSize = figPARAMS.y1FontSize;
                    if figPARAMS.y1labelBold == 1
                        axes1.FontWeight = 'bold';
                    end
                case 'Outside'
                    if (figPARAMS.col(n)==1) || (figPARAMS.row(n)==figPARAMS.subplotSize(1)) % proceed
                        ylabel(figPARAMS.y1label)
                        axes1.FontSize = figPARAMS.y1FontSize;
                        if figPARAMS.y1labelBold == 1
                            axes1.FontWeight = 'bold';
                        end
                    end
                case 'One'
                    if (figPARAMS.row(n)==1) && (figPARAMS.col(n)==1)
                        han = axes(fig,'visible','off'); 
                        han.YLabel.Visible='on';
                        ylabel(han,figPARAMS.y1label);
                        axes1.FontSize = figPARAMS.y1FontSize;
                        if figPARAMS.y1labelBold == 1
                            axes1.FontWeight = 'bold';
                        end
                        set(axes1,'FontWeight','normal')
                    end
            end
        case 'y2'
            switch figPARAMS.labelLoc
                case 'All'
                    setColorbar
                case 'Outside'
                    if (figPARAMS.col(n)==1) || (figPARAMS.row(n)==figPARAMS.subplotSize(1)) % proceed
                        setColorbar
                    end
                case 'One'
                    if (figPARAMS.row(n)==1) && (figPARAMS.col(n)==1)
                        setColorbar
                    end
            end
end