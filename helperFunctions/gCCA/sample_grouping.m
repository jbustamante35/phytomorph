function [H,legendhandle]=sample_grouping(x,y,groupvar,ilab,col,propstr)
% sample_grouping(x,y,groupvar,ilab,col,propstr)
% 
% x,y         Data points to plot
% groupvar    grouping variable, can be char or numerical
% ilab        labels for data points
% col         optional. matrix with colors, length equal to #groups
% propstr     optional. Char array with properties for plot, rows equal to #groups



[N,ii,jj]=unique(groupvar,'rows');

Ngroups=size(N,1);
if Ngroups>64
    Ngroups=64;
end

if length(ilab)==1
    ilab=repmat(ilab,length(x),1);
end

if Ngroups==64 % max number of different colors
   %divide "groupvar" into 64 segments
   if ~isnumeric(N)
       groupvar=jj;
   
   end
   limits=linspace(min(groupvar),max(groupvar),Ngroups);
   groupvar2=[];
   for i=1:length(limits)-1
       idx=find(groupvar>=limits(i) & groupvar<limits(i+1));
       groupvar2(idx)=i;
   end
   idx=find(groupvar>=limits(end));
   groupvar2(idx)=64;
 
   groupvar2=groupvar2';
   N=unique(groupvar2,'rows');
   Ngroups=length(N);
   
else
    groupvar2=groupvar;
end
    

xrange=max(x)-min(x);
H=gcf;
if Ngroups>7
    if ~exist('col','var') | isempty(col)
col=colormap(jet);
    end

colidx=round(linspace(1,size(col,1),Ngroups));
col=col(colidx,:);
else
    if ~exist('col','var')
    col=colormap(lines);
    col=col(1:Ngroups,:);
    end
end

if ~exist('propstr','var') | isempty(propstr)
    propstr=strcat('''color'',[',num2str(col),']');
else
propstr=strcat(propstr,',',strcat('''color'',[',num2str(col),'],''markersize'',12')) ;   
end

hold on;
for i=1:Ngroups
    idx=[];
    for j=1:length(x)
        if isequal(groupvar2(j,:),N(i,:))
            idx=[idx; j];
        end
    end

    h=eval(['plot(x(idx),y(idx),''.'', ' propstr(i,:) ')']);
    
    h=text(x(idx)+xrange/50,y(idx),ilab(idx,:));
    set(h,'color',col(i,:));
end

if size(N,1)<10 & ~isnumeric(N)
legendhandle=legend(num2str(N),'Location','NorthWest');
else
    legendhandle=[];
end

if isnumeric(N) & size(N,1)>9
    H=colorbar;
    
    ticks=get(H,'Ticks');
    ticklabels=linspace(min(groupvar),max(groupvar),length(ticks));
    set(H,'Ticks',ticks);
    set(H,'YTickLabel',num2str(ticklabels'));
end
    
