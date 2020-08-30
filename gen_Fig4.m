%% Generate Figure 4
% Panel a: Bm181021-TDR11 deployment + inset map of animal track
%
% Panel b: Bm190916-TDR14 deployment + inset map of animal track
%
% Each panel shows animal latitude, daily number of feeding lunges
% (circles), and daily number of song calls (triangles). Circles and
% triangles are colored by lunges or calls per hr (night - day).
%
% First section (Lines X-Y) analyzes feeding, vocalizing, and movement 
% behavior from 2018 TDR10 deployment that recorded the transition from 
% foraging to migration (Bm181021-TDR11).
%
% Second section (Lines X-Y) saves a csv to pass to R for statistical 
% analysis of the 2018 individual's behavior before and after the 
% transition to migration.
%
% Third section (Lines X-Y) analyzes feeding, vocalizing, and movement 
% behavior from 2019 TDR10 deployment that recorded the transition from 
% foraging to migration (Bm190916-TDR14).
%
% Fourth section (Lines X-Y) saves a csv to pass to R for statistical 
% analysis of the 2019 individual's behavior before and after the 
% transition to migration.
%
% Fifth section creates main text Figure 4 (panels described above)
%
% Last update: August 30, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; 

% West coast for inset map
thiscoast = shaperead('thiscoast', 'UseGeoCoords', true);

% 2018 TDR10 deployment
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
    lunges_tot(i) = length(lunges);
    [~,call_el] = SolarAzEl(calls+(7/24),36,-121,0);
    callsperhr_day(i) = sum(call_el>0)/(daymins/60);
    callsperhr_night(i) = sum(call_el<=-12)/(nightmins/60);
    callsperhr_dd(i) = (sum(call_el<=0 & call_el>-12))/(ddmins/60); 
    calls_tot(i) = length(calls);
end

% Calculate calls per hour for each solar category, only for days in 
% migratory period

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
disp(['2018 Calls/hr day (migratory) = ',num2str(callsperhr_day_migr)])
callsperhr_night_migr = sum(call_el_migr<=-12)/(nightmins/60);
disp(['2018 Calls/hr night (migratory) = ',num2str(callsperhr_night_migr)])
callsperhr_dd_migr = sum(call_el_migr<=0 & call_el_migr>-12)/(ddmins/60); 
disp(['2018 Calls/hr dusk/dawn (migratory) = ',num2str(callsperhr_dd_migr)])
callsperhr_day_for = sum(call_el_for>0)/(daymins2/60);
disp(['2018 Calls/hr day (foraging) = ',num2str(callsperhr_day_for)])
callsperhr_night_for = sum(call_el_for<=-12)/(nightmins2/60);
disp(['2018 Calls/hr night (foraging) = ',num2str(callsperhr_night_for)])
callsperhr_dd_for = sum(call_el_for<=0 & call_el_for>-12)/(ddmins2/60); 
disp(['2018 Calls/hr dusk/dawn (foraging) = ',num2str(callsperhr_dd_for)])

% vectors for 2018 scatter
x1=days(2:end);
for i=1:length(x1)
    ind1 = find(gps(:,6) >= x1(i),1);
    indend = find(gps(:,6) >= (x1(i)+1),1);
    ind = round(mean([ind1,indend]));
    y1(i) = gps(ind,2) + 0.9;
    y2(i) = gps(ind,2) - 0.9;
end
x1=x1+0.5;
lunges_tot(lunges_tot==0) = NaN;
lungesperhr_night(lunges_tot==0) = NaN;
lungesperhr_day(lunges_tot==0) = NaN;
calls_tot(calls_tot==0) = NaN;
callsperhr_night(calls_tot==0) = NaN;
callsperhr_day(calls_tot==0) = NaN;
sz1=lunges_tot(2:end); 
cl1=(lungesperhr_night-lungesperhr_day)'; cl1=cl1(2:end);
sz2=calls_tot(2:end);
cl2=(callsperhr_night-callsperhr_day)'; cl2=cl2(2:end);

%% 2018 csv for statistics in R
% save csv for use in R
lungesTDR2018 = [lungesperhr_day; lungesperhr_night]'; lungesTDR2018 = lungesTDR2018(2:end,:);
callsTDR2018 = [callsperhr_day; callsperhr_night]'; callsTDR2018 = callsTDR2018(2:end,:);
csvwrite('lungesTDR2018.csv',lungesTDR2018);
csvwrite('callsTDR2018.csv',callsTDR2018);
%% panel b
% 2019 TDR10 deployment
d = 'Bm190916-TDR14';

% Load lunge file and call file
load(['tag_data/TDR10/',d,'/',d,' lunges.mat']);
load(['tag_data/TDR10/',d,'/',d,'_32HzAcc_calls.mat']);

% Subset to only 3-level confidence lunges (highest confidence)
LungeDN = LungeDN(LungeC == 3);

% Get gps fixes and days
gps2 = xlsread(['tag_data/TDR10/',d,'/',d,' GPS.xlsx']);
gpsdv2 = datevec(gps2(:,6));
gpsdays2 = datenum([gpsdv2(:,1),gpsdv2(:,2),gpsdv2(:,3)]);
days2 = starttime:1:datenum('4 Oct 2019');
daysdv2 = datevec(days2);
days2 = datenum([daysdv2(:,1),daysdv2(:,2),daysdv2(:,3)]);

% Calculate daily values for lunges and calls per hour by solar elevation
% category
clear lungesperhr_day lungesperhr_dd lungesperhr_night lunges_tot
clear callsperhr_day callsperhr_dd callsperhr_night calls_tot
for i=1:length(days2)
    idxg = find(gpsdays2==days2(i),1);
    lat = gps2(idxg,2);
    lon = gps2(idxg,3);
    d1 = days2(i);
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
    lunges_tot(i) = length(lunges);
    [~,call_el] = SolarAzEl(calls+(7/24),36,-121,0);
    callsperhr_day(i) = sum(call_el>0)/(daymins/60);
    callsperhr_night(i) = sum(call_el<=-12)/(nightmins/60);
    callsperhr_dd(i) = (sum(call_el<=0 & call_el>-12))/(ddmins/60);
    calls_tot(i) = length(calls);
end

% Calculate calls per hour for each solar category, only for days following
% final lunge (migratory period)

% Behavioral transition date
fdate = datenum('24-Sept-2019 0:00:00');
daymins = 0; nightmins = 0; ddmins = 0;
for i=fdate:1:days2(end)
    d1 = i; d2 = i + 1;
    mins = d1:(1/24/60):(d2-(1/24/60));
    [~,el] = SolarAzEl(mins+(7/24),lat,lon,0);
    daymins = daymins + sum(el>0);
    nightmins = nightmins + sum(el<=-12);
    ddmins = ddmins + 1440 - sum(el>0) - sum(el<=-12);
end
daymins2 = 0; nightmins2 = 0; ddmins2 = 0;
for i=days2(1):1:fdate
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

% Calculate migratory ("migr") and foraging ("for") call rates
% for each solar elevation category. 
callsperhr_day_migr = sum(call_el_migr>0)/(daymins/60);
disp(['2019 Calls/hr day (migratory) = ',num2str(callsperhr_day_migr)])
callsperhr_night_migr = sum(call_el_migr<=-12)/(nightmins/60);
disp(['2019 Calls/hr night (migratory) = ',num2str(callsperhr_night_migr)])
callsperhr_dd_migr = sum(call_el_migr<=0 & call_el_migr>-12)/(ddmins/60); 
disp(['2019 Calls/hr dusk/dawn (migratory) = ',num2str(callsperhr_dd_migr)])
callsperhr_day_for = sum(call_el_for>0)/(daymins2/60);
disp(['2019 Calls/hr day (foraging) = ',num2str(callsperhr_day_for)])
callsperhr_night_for = sum(call_el_for<=-12)/(nightmins2/60);
disp(['2019 Calls/hr night (foraging) = ',num2str(callsperhr_night_for)])
callsperhr_dd_for = sum(call_el_for<=0 & call_el_for>-12)/(ddmins2/60); 
disp(['2019 Calls/hr dusk/dawn (foraging) = ',num2str(callsperhr_dd_for)])

% vectors for 2019 scatter
x2=days2(2:end);
for i=1:length(x2)
    ind1 = find(gps2(:,6) >= x2(i),1);
    indend = find(gps2(:,6) >= (x2(i)+1),1);
    ind = round(mean([ind1,indend]));
    y3(i) = gps2(ind,2) + 0.3;
    y4(i) = gps2(ind,2) - 0.3;
end
x2=x2+0.5;
lunges_tot(lunges_tot==0) = NaN;
lungesperhr_night(lunges_tot==0) = NaN;
lungesperhr_day(lunges_tot==0) = NaN;
calls_tot(calls_tot==0) = NaN;
callsperhr_night(calls_tot==0) = NaN;
callsperhr_day(calls_tot==0) = NaN;
sz3=lunges_tot(2:end);
cl3=(lungesperhr_night-lungesperhr_day)'; cl3=cl3(2:end);
sz4=calls_tot(2:end);
cl4=(callsperhr_night-callsperhr_day)'; cl4=cl4(2:end);

%% 2019 csv for statistics in R
% save csv for use in R
lungesTDR2019 = [lungesperhr_day; lungesperhr_night]'; lungesTDR2019 = lungesTDR2019(2:end,:);
callsTDR2019 = [callsperhr_day; callsperhr_night]'; callsTDR2019 = callsTDR2019(2:end,:);
csvwrite('lungesTDR2019.csv',lungesTDR2019);
csvwrite('callsTDR2019.csv',callsTDR2019);

%% Figure 
% Figure position
figure(1); clf; set(gcf,'position',[200 200 900 600],'color','w');

% Panel positions
P1 = [.1 .57 .85 .42];
P2 = [.1 .07 .85 .42];

% Map position
%Pmap18 = [0.1 .65 .2 .2];

% Font sizes and light gray color
fs = 14; set(0,'DefaultAxesFontSize',fs,'DefaultTextFontsize',fs); 
lg = [.7 .7 .7];

% Diverging colormap
cmap = colormap;
cmap = flipud(cmap);
colormap(cmap);
divmap = NaN(128,3);
divmap(1,:) = [1 0 0];
divmap(64,:) = [1 1 1]; divmap(65,:) = [1 1 1];
divmap(128,:) = [0 0 1];
for i=2:63
    divmap(i,:) = [1 (i-1)/62 (i-1)/62];
end
for i=66:127
    divmap(i,:) = [1-((i-65)/65) 1-((i-65)/65) 1];
end
colormap(divmap)

% Panel a inset map
axes('position',[0.12 0.6 .25 .25]); 
latlim = [27 39.5]; lonlim = [-125 -112];
axesm('MapProjection','Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'PlineLocation',[30 35],'MlineLocation',[-125 -115],'FEdgeColor',lg,...
    'MeridianLabel','on','ParallelLabel','on','MlabelParallel','south',...
    'fontsize',10,'GColor',lg);
geoshow(thiscoast,'FaceColor',lg,'EdgeColor','none');
plotm(gps(:,2),gps(:,3),'k','linewidth',1.5);
framem; gridm; axis off

% Panel a
axes('Position',P1)
plot(gps(:,6),gps(:,2),'k-','Linewidth',2); axis tight; hold on;
h=scatter(x1,y1',sz1,cl1,'Filled','MarkerEdgeColor','k','Linewidth',0.3);
scatter(x1,y2',sz2,cl2,'^','Filled','MarkerEdgeColor','k','Linewidth',0.3)
caxis([-25 25]) 
cb=colorbar;
set(cb,'Ytick',-20:10:20)
ylabel(cb,{'Lunges or calls hr^{-1}';'(night - day)'})

% legend 
hold on
scatter(datenum('9-Nov-2018')+0.5,37,250,'ko')
scatter(datenum('9-Nov-2018')+0.5,38,100,'ko')
scatter(datenum('9-Nov-2018')+0.5,39,10,'ko')
scatter(datenum('18-Nov-2018'),36.8,500,'k^')
scatter(datenum('18-Nov-2018'),38,100,'k^')
scatter(datenum('18-Nov-2018'),38.8,10,'k^')
text(datenum('10-Nov-2018')+0.25,37,'250')
text(datenum('10-Nov-2018')+0.25,38,'100')
text(datenum('10-Nov-2018')+0.25,39,'10')
text(datenum('19-Nov-2018')-0.25,37,'500')
text(datenum('19-Nov-2018')-0.25,38.1,'100')
text(datenum('19-Nov-2018')-0.25,39,'10')
text(datenum('8-Nov-2018')-0.5,40,'Feeding lunges day^{-1}','Fontweight','Bold')
text(datenum('16-Nov-2018')-0.5,40,'Song calls (A+B) day^{-1}','Fontweight','Bold')

datetick('x','mmm dd');  axis tight;
days1 = [days(1)-1;days];
XT = [days1(2:end); max(days1)+1];
ktime = days1(13);
ylim([26 41])
ylabel(['Latitude (',char(176), 'N)']);
yticks([30,35,40])
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
set(gca,'Fontsize',14,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
for T = 1:numel(XT)
    cday = XT(T)+.5; dv = datevec(cday); 
    text(cday,25.5,int2str(dv(3)),'horizontalalignment','center','fontsize',fs); 
end
a = get(gca,'XTickLabel'); 
axis tight;
xlim([datenum('21-Oct-2018'),datenum('23-Nov-2018')]);
frameax
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2);
text(.13,-.1,'October','fontsize',fs); text(.645,-0.1,'November','fontsize',fs)
text(-.1,.98,'A','fontsize',fs,'fontweight','bold')


% Panel b
axes('Position',P2)
plot(gps2(:,6),gps2(:,2),'k-','Linewidth',2); axis tight; hold on;
scatter(x2,y3',sz3,cl3,'Filled','MarkerEdgeColor','k','Linewidth',0.3)
scatter(x2,y4',sz4,cl4,'^','Filled','MarkerEdgeColor','k','Linewidth',0.3)
caxis([-25 25]) 
colormap(divmap)
cb=colorbar;
set(cb,'Ytick',-20:10:20)
ylabel(cb,{'Lunges or calls hr^{-1}';'(night - day)'})

clear days1
datetick('x','mmm dd');  axis tight;
days1 = [days2(1)-1;days2];
XT = [days1(2:end); max(days1)+1];
ktime = days1(10);
ylabel(['Latitude (',char(176), 'N)']);
set(gca,'Fontsize',14,'XTick',XT,'XTickLabel',[],'TickDir','Out','box','off','color','none');
for T = 1:numel(XT)
    cday = XT(T)+.5; dv = datevec(cday); 
    text(cday,31.2,int2str(dv(3)),'horizontalalignment','center','fontsize',fs); 
end
a = get(gca,'XTickLabel'); 
axis tight;
ylim([31.5 37])
ylabel(['Latitude (',char(176), 'N)']);
yticks([33,35,37])
yl = get(gca,'Ylim'); plot(ktime+[0 0],yl,'k--');
xlim([datenum('16-Sep-2019'),datenum('5-Oct-2019')]);
frameax
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2);
text(.18,-.12,'September','fontsize',fs); text(.68,-.12,'October','fontsize',fs)
text(-.1,.98,'B','fontsize',fs,'fontweight','bold')

% Panel b inset map
axes('position',[0.12 0.1 .25 .25]); 
latlim = [27 39.5]; lonlim = [-125 -112];
axesm('MapProjection','Mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'PlineLocation',[30 35],'MlineLocation',[-125 -115],'FEdgeColor',lg,...
    'MeridianLabel','on','ParallelLabel','on','MlabelParallel','south',...
    'fontsize',10,'GColor',lg);
geoshow(thiscoast,'FaceColor',lg,'EdgeColor','none');
plotm(gps2(:,2),gps2(:,3),'k','linewidth',1.5);
framem; gridm; axis off
