function [] = pairWisePlots(X,n,LABEL,curveBank,h,P,toSave,preFix,toWait)
    
    %toSave = '/home/nate/forA/';
    %preFix = 'with_vecs';
    %toWait = true;

    type = {'-','--'};

    if n == 2
        if nargin < 5;close all;else;figure(h);end
        if nargin < 6;P=[];end
        if isempty(P);P = nchoosek((1:size(X,2)),n);end
        CLC = {'r','g','b','k'};
        for e = 1:size(P,1)
            if nargin < 5;figure;else;figure(h);end
            tX = X(:,P(e,:));
            plot(tX(:,1),tX(:,2),'k.');
            xlabel(LABEL{P(e,1)});
            ylabel(LABEL{P(e,2)});
            hold on
            if nargin >= 4
                for bank = 1:numel(curveBank)
                    for c = 1:numel(curveBank{bank})
                        tmpC = curveBank{bank}{c}(:,P(e,:));
                        plot(tmpC(:,1),tmpC(:,2),[CLC{c} type{bank}],'LineWidth',2);
                    end
                end
            end
           
            axis([min(tX(:,1)) max(tX(:,1)) min(tX(:,2)) max(tX(:,2))]);
            
            fileName = [toSave preFix '_2Dplot_' LABEL{P(e,1)} 'VS' LABEL{P(e,2)} '.tif'];
            saveas(gca,fileName);
            
            
            if nargin < 6
                hold off
                if toWait;waitforbuttonpress;end
            else
                hold off
            end
        end
    elseif n == 3
        if nargin < 5;close all;else;figure(h);end
        if nargin < 6;P = nchoosek((1:size(X,2)),n);end
        CLC = {'r','g','b','k'};
        for e = 1:size(P,1)
            if nargin < 5;figure;else;figure(h);end
            tX = X(:,P(e,:));
            plot3(tX(:,1),tX(:,2),tX(:,3),'k.');
            xlabel(LABEL{P(e,1)});
            ylabel(LABEL{P(e,2)});
            zlabel(LABEL{P(e,3)});
            
            hold on
            
            if nargin >= 4
                for bank = 1:numel(curveBank)
                    for c = 1:numel(curveBank{bank})
                        tmpC = curveBank{bank}{c}(:,P(e,:));
                        plot3(tmpC(:,1),tmpC(:,2),tmpC(:,3),[CLC{c} type{bank}],'LineWidth',2);
                    end
                end
            end
            
            fileName = [toSave preFix '_3Dplot_' LABEL{P(e,1)} - LABEL{P(e,2)} '.tif'];
            saveas(gca,fileName);
            
            
            if nargin < 6
                hold off
                 if toWait;waitforbuttonpress;end
                close all
            end
        end
    end
    
    
end