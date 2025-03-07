function updatePairBox(fig,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will update the pairing box to dashed lines for reprocess
% case. 
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(fig); subplot(3,1,3); hold on

ybox = findobj(gca,'Type','line','Color','y');
gbox = findobj(gca,'Type','line','Color','g');
if ~isempty(ybox) 
    if mode == 1 % Dashed for reprocess
        for n = 1:length(ybox)
            set(ybox(n),'LineStyle','--')
        end
    elseif mode == 2 % Gray for und
        for n = 1:length(gbox)
            set(gbox(n),'Color',[179 179 179]./255,'Linewidth',3)
        end
        return
    end
end
if ~isempty(gbox)
    if mode == 2 % Gray for undetermined
        for n = 1:length(gbox)
            set(gbox(n),'Color',[179 179 179]./255,'Linewidth',3)
        end
    end
end
hold off