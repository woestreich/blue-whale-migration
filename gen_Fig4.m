%% Generate Figure 4
% Panel a: TDR10 deployment latitude + inset map of animal track
% Panel b: Daily-binned feeding lunge rate by solar elevation category
% Panel c: Daily-binned night-minus-day call rate
% Last update: May 18, 2020

clear all; close all; 

% West coast for inset map
thiscoast = shaperead('thiscoast', 'UseGeoCoords', true);

% TDR10 deployment
d = 'Bm181021-TDR11';

% Load lunge file and call file
load(['tag_data/TDR10/',d,'/',d,' lunges.mat']);
load(['tag_data/TDR10/',d,'/',d,'_32HzAcc_calls.mat']);

% Subset to only 3-level confidence lunges (highest confidence)
LungeDN = LungeDN(LungeC == 3);

% Get gps fixes and days
gps = xlsread(['tag_data/TDR10/',d,'/',d,' GPS.xlsx']);
gpsdv = datevec(gps(:,6));
gpsdays = datenum([gpsdv(:,1),gpsdv(:,2),gpsdv(:,3)]);
days = starttime:1:datenum('22 Nov 2018');
daysdv = datevec(days);
days = datenum([daysdv(:,1),daysdv(:,2),daysdv(:,3)]);
days = days(2:end-1);

% Calculate daily values for lunges and calls per hour by solar elevation
% category
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
    calls_tot(i) = length(calls);
end

% Calculate calls per hour for each solar category, only for days following
% final lunge (migratory period)

% Behavioral transition date
fdate = datenum('01-Nov-2018 0:00:00');
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

% Calculate and display migratory ("migr") and foraging ("for") call rates
% for each solar elevation category. Reported in Results section of
% manuscript.
callsperhr_day_migr = sum(call_el_migr>0)/(daymins/60);
disp(['Calls/hr day (migratory) = ',callsperhr_day_migr])
callsperhr_night_migr = sum(call_el_migr<=-12)/(nightmins/60);
disp(['Calls/hr night (migratory) = ',callsperhr_night_migr])
callsperhr_dd_migr = sum(call_el_migr<=0 & call_el_migr>-12)/(ddmins/60); 
disp(['Calls/hr dusk/dawn (migratory) = ',callsperhr_dd_migr])
callsperhr_day_for = sum(call_el_for>0)/(daymins2/60);
disp(['Calls/hr day (foraging) = ',callsperhr_day_for])
callsperhr_night_for = sum(call_el_for<=-12)/(nightmins2/60);
disp(['Calls/hr night (foraging) = ',callsperhr_night_for])
callsperhr_dd_for = sum(call_el_for<=0 & call_el_for>-12)/(ddmins2/60); 
disp(['Calls/hr dusk/dawn (foraging) = ',callsperhr_dd_for])

% Figure position
figure(1); clf; set(gcf,'position',[200 200 600 500],'color','w');

% Panel positions
px = .125; pw = .85;
ph = .285; py = .96-[1:3]*ph-[0:2]*.01;
P1 = [px py(1) pw ph]; P2 = [px py(2) pw ph]; P3 = [px py(3) pw ph];

% Map position
Pmap18 = [0.1 .65 .2 .2];

% Font sizes
fs = 12; set(0,'DefaultAxesFontSize',fs,'DefaultTextFontsize',fs); 

% Grays 
lg = [.7 .7 .7]; dg = [.45 .45 .45];

% X ticks
days = [days(1)-1;days]+0.5; 
XT = [days; max(days)+1]-.5;
ktime = days(11)+0.5;

% Panel a
axes('position',P1);
plot(gps(:,6),gps(:,2),'k-','Linewidth',2); axis tight;
ylabel(['Latitude (',char(176), 'N)']);
ylim([27.5 40]); xlim([datenum('21-Oct-2018'),datenum('23-Nov-2018')]);
hold on
yticks([30,35,40])
set(gca,'Fontsize',12,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
frameax
text(-0.1,0.96,'a','Fontsize',12,'fontweight','bold')

% Panel b
axes('position',P2);
lunges = [lungesperhr_night; lungesperhr_dd; lungesperhr_day]';
lunges = [NaN, NaN, NaN;lunges];
b1=bar(days,lunges,'hist'); 
datetick('x','mmm dd'); axis tight;
set(b1(1),'FaceColor',[0.4 0.4 0.4])
set(b1(1),'EdgeColor','k')
set(b1(2),'FaceColor',[0.8,0.8,0.8])
set(b1(2),'EdgeColor','k')
set(b1(3),'FaceColor','w')
set(b1(3),'EdgeColor','k')
hold on;
ylabel({'Feeding rate';'(lunges hr^{-1})'});
ylim([0 22]); xlim([datenum('21-Oct-2018'),datenum('23-Nov-2018')]);
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
legend('Night','Dusk/Dawn','Day')
set(gca,'Fontsize',12,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
frameax
text(-0.1,0.96,'b','Fontsize',12,'fontweight','bold')

% Panel c
axes('position',P3);
b2=bar(days,[NaN,callsperhr_night - callsperhr_day],0.8);
datetick('x','mmm dd');  axis tight;
co = get(gca,'colororder'); mb = co(1,:);
set(b2(1),'FaceColor',mb)
set(b2(1),'EdgeColor',mb)
hold on;
% Annotation
lx = datenum([2018 11 22 12 0 0]); ly = [1 22]; 
dg = [0 0 0];
plot(lx+[0 0],ly,'k','linewidth',2,'color',dg); 
plot(lx,max(ly),'k^','markersize',7,'color',dg,'markerfacecolor',dg);
text(lx-4.25,mean(ly)+1,'More song','color',dg);
text(lx-4.25,mean(ly)-4,'during night','color',dg);
plot(lx+[0 0],-ly,'k','linewidth',2,'color',dg); 
plot(lx,-max(ly),'kv','markersize',7,'color',dg,'markerfacecolor',dg);
text(lx-4.25,-mean(ly)+1,'More song','color',dg);
text(lx-4.25,-mean(ly)-3.5,'during day','color',dg);
set(gca,'Fontsize',12,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
for T = 1:numel(XT);
    cday = XT(T)+.5; dv = datevec(cday); 
    text(cday,-31,int2str(dv(3)),'horizontalalignment','center','fontsize',fs); 
end
a = get(gca,'XTickLabel');  
ylabel({'Calls night - day';'(calls hr^{-1})'});
ylim(28*[-1 1]) 
xlim([datenum('21-Oct-2018'),datenum('23-Nov-2018')]);
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
frameax
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2);
text(-0.1,0.96,'c','Fontsize',12,'fontweight','bold')

hold on;

% Place map on Panel a
axes('position',[px+.05 py(1)+.025 .2 .2]); 
latlim = [27 39.5]; lonlim = [-125 -112];
axesm('MapProjection','Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'PlineLocation',[30 35],'MlineLocation',[-125 -115],'FEdgeColor',lg,...
    'MeridianLabel','on','ParallelLabel','on','MlabelParallel','south',...
    'fontsize',10,'GColor',lg);
geoshow(thiscoast,'FaceColor',lg,'EdgeColor','none');
plotm(gps(:,2),gps(:,3),'k','linewidth',1.5);
framem; gridm; axis off

% Axis positioning and month labels on x axis
axes('position',[0 0 1 1]); set(gca,'visible','off');
text(.23,.04,'October','fontsize',12); text(.645,.04,'November','fontsize',12)
