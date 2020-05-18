%% Frame axes in matlab plot
%  Written by John Ryan; MBARI
g=get(gca,'position'); axes('position',g); 
set(gca,'color','none','box','on','XTick',[],'YTick',[]); %,'linewidth',1);

