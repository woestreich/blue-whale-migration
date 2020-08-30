%% Generate Figure 3
% Panel a: Image of CATS tag deployment on blue whale
%
% Panel b: Temporal coverage of CATS and TDR10 deployments analyzed in this
% study
%
% Panel c: Map of CATS tag deployment and MARS hydrophone locations  
%
% Panel d: Example dive profile and ethogram from CATS deployment
%
% Panel e: Aggregated 2017-2019 results for call rates
%
% Panel f: Aggregated 2017-2019 results for feeding rates
% 
% Last update: August 30, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; 

% aggregate tag info for Panel b
% aggregate CATS 2017 info
catsdir = dir('tag_data/CATS/'); catsdir17 = catsdir(4:7);
for i=1:length(catsdir17)
    load(['tag_data/CATS/',catsdir17(i).name,'/',catsdir17(i).name,'_taginfo.mat'])
    dn = DN(tagon);
    tags17(i,1) = dn(1);
    tags17(i,2) = dn(end);
end
tags17 = tags17 - datenum('01-Jan-2017 00:00:00');
yr17 = [2016.7, 2016.7; 2016.9, 2016.9; 2017.1, 2017.1; 2017.3, 2017.3];

% aggregate CATS 2018 info
catsdir = dir('tag_data/CATS/'); catsdir18 = catsdir(8:13);
for i=1:length(catsdir18)
    load(['tag_data/CATS/',catsdir18(i).name,'/',catsdir18(i).name,'_taginfo.mat'])
    dn = DN(tagon);
    tags18(i,1) = dn(1);
    tags18(i,2) = dn(end);
end
tags18 = tags18 - datenum('01-Jan-2018 00:00:00');
yr18 = [2017.6, 2017.6; 2017.76, 2017.76; 2017.92, 2017.92; 2018.08, 2018.08; 2018.24, 2018.24; 2018.4, 2018.4];

% aggregate CATS 2019 info
catsdir = dir('tag_data/CATS/'); catsdir19 = catsdir(14:16);
for i=1:length(catsdir19)
    load(['tag_data/CATS/',catsdir19(i).name,'/',catsdir19(i).name,'_taginfo.mat'])
    dn = DN(tagon);
    tags19(i,1) = dn(1);
    tags19(i,2) = dn(end);
end
tags19 = tags19 - datenum('01-Jan-2019 00:00:00');
yr19 = [2018.8, 2018.8; 2019, 2019; 2019.2, 2019.2];

% TDR10 2018
load('tag_data/TDR10/Bm181021-TDR11/Bm181021-TDR11_taginfo.mat')
dn = DN(tagon);
tdr18 = [dn(1) dn(end)];
tdr18 = tdr18 - datenum('01-Jan-2018 00:00:00');

% TDR10 2019
load('tag_data/TDR10/Bm190916-TDR14/Bm190916-TDR14_taginfo.mat')
dn = DN(tagon);
tdr19 = [dn(1) dn(end)];
tdr19 = tdr19 - datenum('01-Jan-2019 00:00:00');

clear DN
% Blue whale CATS tag deployment image (for panel a)
bi = imread('tag_data/blue_zoom.jpg'); 

% Set up figure positioning
figure(1); clf; set(gcf,'position',[200 200 700 700],'color','w');
fs = 11;
px = .09; pw = .89; pw2 = .44*pw;
Apos = [px .73 pw .31];
Bpos = [px .6 pw-0.11 .1];
Cpos = [px+pw-0.105 .6 .12 .12];
Dpos = [px .3 pw .24];
Epos = [px .03 pw2 .24];
Fpos = [px+pw2+.107 .03 pw2 .24];

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
text(-250,150,'A','fontsize',12,'fontweight','bold')

% Panel b
axes('position',Bpos);
figure(1)
plot(tags17',yr17','k-o','MarkerFaceColor','k','MarkerSize',4); hold on;
plot(tags18',yr18','k-o','MarkerFaceColor','k','MarkerSize',4);
plot(tdr18,[2018,2018],'r-o','MarkerFaceColor','r','MarkerSize',4);
plot(tags19',yr19','k-o','MarkerFaceColor','k','MarkerSize',4);
plot(tdr19,[2019,2019],'r-o','MarkerFaceColor','r','MarkerSize',4);
set(gca,'Xtick',[212,243,273,304,334],'Xlim',[212 334])
set(gca,'XtickLabel',{'                                          Aug',...
    '                                          Sep',...
    '                                          Oct',...
    '                                          Nov',''},'Xaxislocation','top','Fontsize',12)
set(gca,'Ytick',[],'Ylim',[2016.5 2019.5])
text(199.5,2019.8,'B','fontsize',12,'fontweight','bold');
text(334.5,2019.8,'C','fontsize',12,'fontweight','bold');
plot([212,334],[2017.5,2017.5],'k--')
plot([232,334],[2018.5,2018.5],'k--')
text(205.8,2017,'2017','fontsize',11); 
text(205.8,2018,'2018','fontsize',11);
text(205.8,2019,'2019','fontsize',11);
patch([212,232,232,212],[2017.9,2017.9,2019.1,2019.1],'w')
plot([221,231],[2018.8,2018.8],'k-o','Markerfacecolor','k','MarkerSize',4);
plot([221,231],[2018.2,2018.2],'r-o','Markerfacecolor','r','MarkerSize',4);
text(212.7,2018.8,'CATS','fontsize',11)
text(212.7,2018.2,'TDR10','fontsize',11)
box on

% Panel c
cats_locs(:,1) = [36.8070; 36.8072; 36.8015; 36.8027; 36.8987; 36.9037; ...
    36.7827; 36.8117; 36.8025; 36.7992; 36.6715; 36.5050; 36.7146];
cats_locs(:,2) = [-122.0379; -122.0817; -122.0750; -122.0730; -122.2640; ...
    -122.2682; -122.0543; -122.1383; -122.1308; -122.0912; -122.0081; ...
    -121.9694; -121.9818];
AX = [-122.4 -121.7 36.45 37.07];
axes('position',Cpos);
axesm('MapProjection','Mercator','MapLatLimit',AX([3 4]),'MapLonLimit',AX([1 2]),...
    'fontsize',12,'GColor','k','FEdgeColor',[.5 .5 .5],'FontColor','k','FLineWidth',1);
framem; gridm; axis off
% West coast for map
thiscoast = shaperead('thiscoast', 'UseGeoCoords', true);
geoshow(thiscoast,'FaceColor',[.7 .7 .7],'EdgeColor','k');
plotm(cats_locs(:,1),cats_locs(:,2),'k.'); hold on;
% MARS hydrophone location
load acoustic_data/MARS_hydrophone_location; MARS.hlat = hloc.lat; MARS.hlon = hloc.lon;
plotm(MARS.hlat,MARS.hlon,'ko'); textm(MARS.hlat-.07,MARS.hlon-.2,'MARS','fontsize',10); 

% Panel d
% Data for Figure 3d
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
axes('position',Dpos);
p1=patch([x(idxn),x(idxn),DN(end),DN(end)],[[-100,15],[15,-100]],[0.4,0.4,0.4],'HandleVisibility','off');
hold on
p2=patch([x(idxdd),x(idxdd),x(idxn),x(idxn)],[[-100,15],[15,-100]],[0.8,0.8,0.8],'HandleVisibility','off');
green = [0 0.8 0.2];
cyan = [0.2 0.6 0.9];
plot(DN,p*(-1),'color','k'); hold on;...
    plot(AB_times_DN,AB_depths*-1,'.','color',green,'MarkerSize',10);...
    plot(DN(LungeI),p(LungeI)*(-1),'.','color',cyan,'MarkerSize',10)
set(gca,'Fontsize',fs)
[h,icons] = legend({'Depth (m)','Song call','Feeding lunge'},'Position',[0.72,0.33,0.18,0.13],'Fontsize',11); 
icons = findobj(icons,'type','line');
set(icons,'MarkerSize',20)
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
text(DN(1)-0.08,20,'D','Fontsize',12,'fontweight','bold')
axis tight
box on
ylim([-99,15])
text(datenum('Aug 14 2017 13:30'),-92,'Day','Fontsize',fs-1)
text(datenum('Aug 14 2017 20:00'),-92,'Dusk','Fontsize',fs-1)
text(datenum('Aug 15 2017 0:40'),-92,'Night','Fontsize',fs-1,'Color','w')

% Load in aggregated tag data for panels e and f
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

% Panel e
axes('position',Epos);
p1=bar(bar_data_a);
set(p1(1),'FaceColor','w'); set(p1(2),'FaceColor',[0.4 0.4 0.4]); set(p1(3),'FaceColor',[0.8,0.8,0.8]);
set(gca, 'XTick', [1 2 3], 'XLim', [.5 3.5], 'YLim', [0 12])
set(gca, 'XTickLabel', {'2017' '2018' '2019'})
set(gca, 'Fontsize', fs)
ylabel({'Song call rate' ; '(A + B calls hour^{-1})'})
text(-0.1,12.5,'E','Fontsize',12,'fontweight','bold')
% Custom legend to fit in a smaller space
lx = 2.65; ly = 11.1; dy = 1.25;
hold on;
plot(lx,ly,'ks','markersize',10); 
text(lx+.1,ly,'Day','fontsize',fs);
plot(lx,ly-dy,'ks','markersize',10,'markerfacecolor',[.4 .4 .4]); 
text(lx+.1,ly-dy,'Night','fontsize',fs);
plot(lx,ly-2*dy,'ks','markersize',10,'markerfacecolor',[.8 .8 .8]);
text(lx+.1,ly-2*dy,'Dusk/Dawn','fontsize',fs);

% Panel f
axes('position',Fpos);
p2=bar(bar_data_b);
set(p2(1),'FaceColor','w'); set(p2(2),'FaceColor',[.4 .4 .4]); set(p2(3),'FaceColor',[0.8,0.8,0.8]);
set(gca, 'XTick', [1 2 3], 'XLim', [0.5 3.5], 'YLim', [0 22])
set(gca, 'XTickLabel', {'2017' '2018' '2019'})
set(gca, 'Fontsize', fs)
ylabel({'Feeding rate' ; '(lunges hour^{-1})'})
text(-0.1,23,'F','Fontsize',12,'fontweight','bold')


