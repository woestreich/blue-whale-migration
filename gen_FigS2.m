% Extended data figure 2

clear all; close all; 

% load in and process TDR10 deployment data
d = 'Bm190916-TDR14';

%load lunge file and call file
load(['tag_data/TDR10/',d,'/',d,' lunges.mat']);
load(['tag_data/TDR10/',d,'/',d,'_32HzAcc_calls.mat']);

%subset to only 3-level confidence lunges (for sure lunges)
LungeDN = LungeDN(LungeC == 3);

% get gps fixes and days
gps = xlsread(['tag_data/TDR10/',d,'/',d,' GPS.xlsx']);
gpsdv = datevec(gps(:,6));
gpsdays = datenum([gpsdv(:,1),gpsdv(:,2),gpsdv(:,3)]);
days = starttime:1:datenum('4 Oct 2019');
daysdv = datevec(days);
days = datenum([daysdv(:,1),daysdv(:,2),daysdv(:,3)]);
days = days(2:end-1);

for i=1:length(days)
    idxg = find(gpsdays==days(i),1);
    lat = gps(idxg,2);
    lon = gps(idxg,3);
    d1 = days(i);
    d2 = d1 + 1;
    idxl = LungeDN >= d1 & LungeDN < d2;
    lunges = LungeDN(idxl);
    idxc = callDN >= d1 & callDN < d2;
    calls = callDN(idxc);
    mins = d1:(1/24/60):(d2-(1/24/60));
    [~,el] = SolarAzEl(mins+(7/24),lat,lon,0);
    daymins = sum(el>0);
    nightmins = sum(el<=-12);
    ddmins = 1440 - daymins - nightmins;
    [~,lunge_el] = SolarAzEl(lunges+(7/24),36,-121,0);
    lungesperhr_day(i) = sum(lunge_el>0)/(daymins/60);
    lungesperhr_night(i) = sum(lunge_el<=-12)/(nightmins/60);
    lungesperhr_dd(i) = (sum(lunge_el<=0 & lunge_el>-12))/(ddmins/60);
    [~,call_el] = SolarAzEl(calls+(7/24),36,-121,0);
    callsperhr_day(i) = sum(call_el>0)/(daymins/60);
    callsperhr_night(i) = sum(call_el<=-12)/(nightmins/60);
    callsperhr_dd(i) = (sum(call_el<=0 & call_el>-12))/(ddmins/60);
    
end

%calculate calls per hour for each solar category, only for days following
%final lunge (migratory period)
fdate = datenum('24-Sept-2019 0:00:00');
daymins = 0; nightmins = 0; ddmins = 0;
for i=fdate:1:days(end)
    d1 = i; d2 = i + 1;
    mins = d1:(1/24/60):(d2-(1/24/60));
    [~,el] = SolarAzEl(mins+(7/24),lat,lon,0);
    daymins = daymins + sum(el>0);
    nightmins = nightmins + sum(el<=-12);
    ddmins = ddmins + 1440 - sum(el>0) - sum(el<=-12);
end
daymins2 = 0; nightmins2 = 0; ddmins2 = 0;
for i=days(1):1:fdate
    d1 = i; d2 = i + 1;
    mins = d1:(1/24/60):(d2-(1/24/60));
    [~,el] = SolarAzEl(mins+(7/24),lat,lon,0);
    daymins2 = daymins2 + sum(el>0);
    nightmins2 = nightmins2 + sum(el<=-12);
    ddmins2 = ddmins2 + 1440 - sum(el>0) - sum(el<=-12);
end
calls_migration = callDN(callDN >= fdate);
calls_forage = callDN(callDN < fdate);
[~,call_el_migr] = SolarAzEl(calls_migration+(7/24),36,-121,0);
[~,call_el_for] = SolarAzEl(calls_forage+(7/24),36,-121,0);
callsperhr_day_migr = sum(call_el_migr>0)/(daymins/60)
callsperhr_night_migr = sum(call_el_migr<=-12)/(nightmins/60)
callsperhr_dd_migr = sum(call_el_migr<=0 & call_el_migr>-12)/(ddmins/60) 
callsperhr_day_for = sum(call_el_for>0)/(daymins2/60)
callsperhr_night_for = sum(call_el_for<=-12)/(nightmins2/60)
callsperhr_dd_for = sum(call_el_for<=0 & call_el_for>-12)/(ddmins2/60)

thiscoast = shaperead('thiscoast', 'UseGeoCoords', true);

figure(1); clf; set(gcf,'position',[200 200 600 500],'color','w');

%Panel positions
px = .125; pw = .85;
ph = .285; py = .96-[1:3]*ph-[0:2]*.01;
P1 = [px py(1) pw ph]; P2 = [px py(2) pw ph]; P3 = [px py(3) pw ph];
%Map position
Pmap18 = [0.1 .65 .2 .2];

fs = 12; set(0,'DefaultAxesFontSize',fs,'DefaultTextFontsize',fs); 

% grays 
lg = [.7 .7 .7]; dg = [.45 .45 .45];

% X ticks
days = [days(1)-1;days]+0.5; 
XT = [days; max(days)+1]-.5;
ktime = days(8)+0.5;

%Panel A
axes('position',P1);
plot(gps(:,6),gps(:,2),'k-','Linewidth',2); axis tight;
ylabel(['Latitude (',char(176), 'N)']);
ylim([32 37]); xlim([datenum('16-Sep-2019'),datenum('5-Oct-2019')]);
hold on
yticks([33,35,37])
set(gca,'Fontsize',12,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
frameax
% title({'Blue whale - tag 11';'2018'},'Fontsize',14)
text(-0.1,0.96,'a','Fontsize',12,'fontweight','bold')

%Panel B
axes('position',P2);
lunges = [lungesperhr_night; lungesperhr_dd; lungesperhr_day]';
lunges = [NaN, NaN, NaN;lunges];
b1=bar(days,lunges,'hist'); 
% b1=bar(days,lunges,'stacked'); 
datetick('x','mmm dd'); axis tight;
set(b1(1),'FaceColor',[0.4 0.4 0.4])
set(b1(1),'EdgeColor','k')
set(b1(2),'FaceColor',[0.8,0.8,0.8])
set(b1(2),'EdgeColor','k')
set(b1(3),'FaceColor','w')
set(b1(3),'EdgeColor','k')
hold on;
% xline(days(10)+0.5,'k--','Linewidth',2);
ylabel({'Feeding rate';'(lunges hr^{-1})'});
ylim([0 22]); xlim([datenum('16-Sep-2019'),datenum('5-Oct-2019')]); 
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
set(gca,'Xtick',days(1)-0.5:1:days(end)+0.5)
legend('Night','Dusk/Dawn','Day')
set(gca,'Fontsize',12,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
frameax
text(-0.1,0.96,'b','Fontsize',12,'fontweight','bold')

%Panel C
axes('position',P3);
b2=bar(days,[NaN,callsperhr_night - callsperhr_day],0.8);
datetick('x','mmm dd');  axis tight;
co = get(gca,'colororder'); mb = co(1,:);
set(b2(1),'FaceColor',mb)
set(b2(1),'EdgeColor',mb)
hold on;
% Annotation, old school
lx = datenum([2019 10 4 18 0 0]); ly = [1 22]; 
dg = [0 0 0];
plot(lx+[0 0],ly,'k','linewidth',2,'color',dg); 
plot(lx,max(ly),'k^','markersize',7,'color',dg,'markerfacecolor',dg);
text(lx-2.5,mean(ly)+3,'More song','color',dg);
text(lx-2.5,mean(ly)-2.5,'during night','color',dg);
plot(lx+[0 0],-ly,'k','linewidth',2,'color',dg); 
plot(lx,-max(ly),'kv','markersize',7,'color',dg,'markerfacecolor',dg);
text(lx-2.5,-mean(ly)-2.8,'More song','color',dg);
text(lx-2.5,-mean(ly)-7,'during day','color',dg);


set(gca,'Fontsize',12,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
for T = 1:numel(XT);
    cday = XT(T)+.5; dv = datevec(cday); 
    text(cday,-31,int2str(dv(3)),'horizontalalignment','center','fontsize',fs); 
end

a = get(gca,'XTickLabel');  
ylabel({'Calls night - day';'(calls hr^{-1})'});
ylim(28*[-1 1]) 
xlim([datenum('16-Sep-2019'),datenum('5-Oct-2019')]);
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
frameax
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2);
text(-0.1,0.96,'c','Fontsize',12,'fontweight','bold')

%annotation arrows
hold on;
% annotation('arrow',[0.88,0.88],[0.24,0.33])
% annotation('arrow',[0.88,0.88],[0.23,0.14])
% text(0.88,0.75,{'More song';'in night'},'Fontsize',10)
% text(0.88,0.2,{'More song';'in day'},'Fontsize',10)

%place map on Panel A
% axes('position',Pmap18); imshow('map2018.png')
axes('position',[px+.675 py(1)+.075 .2 .2]); 
latlim = [27 39.5]; lonlim = [-125 -112];
axesm('MapProjection','Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'PlineLocation',[30 35],'MlineLocation',[-125 -115],'FEdgeColor',lg,...
    'MeridianLabel','on','ParallelLabel','on','MlabelParallel','south',...
    'fontsize',10,'GColor',lg);
geoshow(thiscoast,'FaceColor',lg,'EdgeColor','none');
plotm(gps(:,2),gps(:,3),'k','linewidth',1.5);
framem; gridm; axis off

axes('position',[0 0 1 1]); set(gca,'visible','off');
text(.41,.04,'September','fontsize',12); text(.85,.04,'October','fontsize',12)

