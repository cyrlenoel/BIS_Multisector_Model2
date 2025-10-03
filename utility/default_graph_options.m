% default_graph_options
%
% Selection of default graph options
% 
% =========================================================================


box off ; 
ax          = gca ; 
ax.YGrid    = 'on' ;
set(get(gcf,'CurrentAxes'),'FontSize',15, 'LineWidth',1.5)    


if exist('graph_options') ==1 
    
    if isfield(graph_options,'xlim') 
        ax.XLim = graph_options.xlim ;
    end
    
    if isfield(graph_options,'xticks') 
        ax.XTick = graph_options.xticks ;
    end
    
    if isfield(graph_options,'xticklabels') 
        ax.XTickLabel = graph_options.xticklabels ;
    end
       
    if isfield(graph_options,'ylabel') 
        ax.YLabel.String = cellstr(graph_options.ylabel) ;
    end

    if isfield(graph_options,'xlabel') 
        ax.XLabel.String = cellstr(graph_options.xlabel) ;
    end    
    
    if isfield(graph_options,'zeroline') 
        line([ax.XLim(1) ax.XLim(end)], [0 0], 'LineWidth',1, 'Color', 'k') ;
    end
    
    if isfield(graph_options,'ylim') 
        ax.YLim = graph_options.ylim ;
    end
    
    if isfield(graph_options,'yticks') 
        ax.YTick = graph_options.yticks ;
    end

    if isfield(graph_options,'title') 
        ax.Title.String = cellstr(graph_options.title) ;
    end
    
    
end
    
    

