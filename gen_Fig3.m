%% Generate Figure 3
% Panel a: Image of CATS tag deployment on blue whale
% Panel b: Example dive profile and ethogram from CATS deployment
% Panel c: Aggregated 2017-2019 results for call rates
% Panel d: Aggregated 2017-2019 results for feeding rates
% Last update: May 18, 2020

clear all; close all; 

% Blue whale CATS tag deployment image (for panel a)
bi = imread('tag_data/blue_zoom.jpg'); 

% Data for Figure 3b
deployment = 'bw170814-51';
load(['tag_data/CATS/',deployment,'/',deployment,'_taginfo.mat'])

% load in corresponding lunge detection file
load(['tag_data/CATS/',deployment,'/',deployment,'lunges.mat'])

% load in phrase and unit times and types
phrasetime_fname = ['tag_data/CATS/',deployment,'/phrasetimes.csv'];
phrasetime = csvread(phrasetime_fname, 1, 0);

calltime_fname = ['tag_data/CATS/',deployment,'/calltimes.csv'];
calltime = csvread(calltime_fname, 1, 0);

phrasetype_fname = ['tag_data/CATS/',deployment,'/phrasetypes.csv'];
delimiter = '';
formatSpec = '%s%[^\n\r]';
fileID = fopen(phrasetype_fname,'r');
phrase = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
phrasetype = phrase{1};
phrasetype = phrasetype(2:length(phrasetype));
fclose(fileID);

calltype_fname = ['tag_data/CATS/',deployment,'/callunits.csv'];
callunit = table2array(readtable(calltype_fname)); 

% sync up phrase time DN's with PRH file DN's (so depth value can be assigned to each call
calltime=DN(1)+((1/86400).*calltime);
for i=1:length(calltime)
    tdiffs=calltime(i) - DN;
    [mindiff, mindiff_ind]=min(abs(tdiffs));
    calltime_DN(i)=DN(mindiff_ind);
    calldepths(i)=p(mindiff_ind);
end

phrasetime=DN(1)+((1/86400).*phrasetime);
for i=1:length(phrasetime)
    tdiffs=phrasetime(i) - DN;
    [mindiff, mindiff_ind]=min(abs(tdiffs));
    phrasetime_DN(i)=DN(mindiff_ind);
    phrasedepths(i)=p(mindiff_ind);
end
A_times_DN = calltime_DN(strcmp(callunit,'A'));
B_times_DN = calltime_DN(strcmp(callunit,'B'));
AB_times_DN = [A_times_DN, B_times_DN];
A_depths = calldepths(strcmp(callunit,'A'));
B_depths = calldepths(strcmp(callunit,'B'));
AB_depths = [A_depths, B_depths];

% find day, night, dusk/dawn indices for this deployment
x = DN(1):(1/24/60):(DN(end)-(1/24/60));
for i=1:length(x)
    [~,xel(i)] = SolarAzEl(x(i)+(7/24),36.789,-122.078,0);
end 
idxn = find(xel<-12,1);
idxdd = find(xel<0,1);

% Set up figure positioning
figure(1); clf; set(gcf,'position',[200 200 600 450],'color','w');
fs = 11;
px = .09; pw = .89; pw2 = .44*pw;
Apos = [px .67 pw .31];
Bpos = [px .33 pw .25];
Cpos = [px .05 pw2 .25];
Dpos = [px+pw2+.107 .05 pw2 .25];

% Panel a, blue whale image
% Avoid distorting the aspect ratio by relating the figure and image aspect
% ratios
axes('position',Apos);
image(bi); set(gca,'Xtick',[],'Ytick',[],'box','on'); 
hold on;
lx = 1130; ly = 560;
plot(lx,ly,'w^','markersize',8,'markerfacecolor','w');
plot(lx+[0 0],ly+[0 100],'w','linewidth',2);
text(lx,ly+140,'CATS tag','fontsize',12,'fontweight','bold','color','w','Horizontalalignment','center');
text(-250,25,'a','fontsize',12,'fontweight','bold');

% Panel b
axes('position',Bpos);
p1=patch([x(idxn),x(idxn),DN(end),DN(end)],[[-100,15],[15,-100]],[0.4,0.4,0.4],'HandleVisibility','off');
hold on
p2=patch([x(idxdd),x(idxdd),x(idxn),x(idxn)],[[-100,15],[15,-100]],[0.8,0.8,0.8],'HandleVisibility','off');
green = [0 0.8 0.2];
cyan = [0.2 0.6 0.9];
plot(DN,p*(-1),'color','k'); hold on;...
    plot(AB_times_DN,AB_depths*-1,'.','color',green,'MarkerSize',8);...
    plot(DN(LungeI),p(LungeI)*(-1),'.','color',cyan,'MarkerSize',8)
set(gca,'Fontsize',fs)
legend({'Depth (m)','Song call','Feeding lunge'},'Position',[0.72,0.365,0.18,0.13],'Fontsize',11); 
yl = get(gca,'ylim');
set(gca,'xticklabel',datestr((get(gca,'xtick').'),'HH:MM'),'Xaxislocation','top')
ytks = [-100:25:0];
yticks(ytks)
set(gca,'yticklabel',abs(ytks))
xtks = [datenum('Aug 14 2017 8:00:00'):1/6:datenum('Aug 15 2017 4:00:00')];
xticks(xtks)
set(gca, 'Layer', 'top')
datetick('x','mm/dd HH:MM','keeplimits','keepticks')
ylabel('Depth (m)')
xlabel('Local time')
text(DN(1)-0.08,10,'b','Fontsize',12,'fontweight','bold')
axis tight
box on
ylim([-99,15])
text(datenum('Aug 14 2017 13:30'),-92,'Day','Fontsize',fs-1)
text(datenum('Aug 14 2017 20:00'),-92,'Dusk','Fontsize',fs-1)
text(datenum('Aug 15 2017 0:40'),-92,'Night','Fontsize',fs-1,'Color','w')

% Load in aggregated tag data for panels c and d
bar_data_a = ones(3,3);
bar_data_b = ones(3,3);
load 'tag_data/CATS_TDR_diel_2017.mat'
bar_data_a(1,:) = callrate;
bar_data_b(1,:) = lungerate;
load 'tag_data/CATS_TDR_diel_2018.mat'
bar_data_a(2,:) = callrate;
bar_data_b(2,:) = lungerate;
load 'tag_data/CATS_TDR_diel_2019.mat'
bar_data_a(3,:) = callrate;
bar_data_b(3,:) = lungerate;

% Panel c
axes('position',Cpos);
p1=bar(bar_data_a);
set(p1(1),'FaceColor','w'); set(p1(2),'FaceColor',[0.4 0.4 0.4]); set(p1(3),'FaceColor',[0.8,0.8,0.8]);
set(gca, 'XTick', [1 2 3], 'XLim', [.5 3.5], 'YLim', [0 12])
set(gca, 'XTickLabel', {'2017' '2018' '2019'})
set(gca, 'Fontsize', fs)
ylabel({'Song call rate' ; '(A + B calls hour^{-1})'})
text(-0.1,12.5,'c','Fontsize',12,'fontweight','bold')
% Custom legend to fit in a smaller space
lx = 2.65; ly = 11.1; dy = 1.25;
hold on;
plot(lx,ly,'ks','markersize',10); 
text(lx+.1,ly,'Day','fontsize',fs);
plot(lx,ly-dy,'ks','markersize',10,'markerfacecolor',[.4 .4 .4]); 
text(lx+.1,ly-dy,'Night','fontsize',fs);
plot(lx,ly-2*dy,'ks','markersize',10,'markerfacecolor',[.8 .8 .8]);
text(lx+.1,ly-2*dy,'Dusk/Dawn','fontsize',fs);

% Panel d
axes('position',Dpos);
p2=bar(bar_data_b);
set(p2(1),'FaceColor','w'); set(p2(2),'FaceColor',[.4 .4 .4]); set(p2(3),'FaceColor',[0.8,0.8,0.8]);
set(gca, 'XTick', [1 2 3], 'XLim', [0.5 3.5], 'YLim', [0 22])
set(gca, 'XTickLabel', {'2017' '2018' '2019'})
set(gca, 'Fontsize', fs)
ylabel({'Lunge rate' ; '(lunges hour^{-1})'})
text(-0.1,23,'d','Fontsize',12,'fontweight','bold')