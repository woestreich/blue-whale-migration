%% Generate Figure 2
% Panel a: modeled received level results for MARS hydrophone
% Panel b: average annual cycle of call index (CI)
% Panel c: average annual cycle of delta CI
% Last update: May 18, 2020

clear all; close all; 

% Daily power spectral density binned by solar elevation 
load 'acoustic_data/MARS_SpectrumLevel_Daily_bySolarElevation.mat'
    
% Handle for month
dv = datevec(D.time); mo = dv(:,2); 

% Climatology across all solar elevation bins
sm = sum(D.sm,3); ct = sum(D.ct,3);
for M = 1:12;
    idx = find(mo == M);  
    L.ltsa(:,M) = sum(sm(:,idx),2) ./ sum(ct(:,idx),2);
end
L.freq = D.freq; L.time = [1:12]; 
CI.all = call_index(L);
    
% Climatology for each solar elevation bin
sn = {'n','dd','d'};
for S = 1:3;
    sm = D.sm(:,:,S); ct = D.ct(:,:,S);
    clear L;
    for M = 1:12;
        idx = find(mo == M);
        L.ltsa(:,M) = sum(sm(:,idx),2) ./ sum(ct(:,idx),2);
    end
    L.freq = D.freq; L.time = [1:12];
    eval(['CI.' sn{S} ' = call_index(L);']);
end

% Received level model output
RL = load('acoustic_data/RL_Sc01_November.mat');

% West coast for map
thiscoast = shaperead('thiscoast', 'UseGeoCoords', true);

% MARS hydrophone location
load acoustic_data/MARS_hydrophone_location; MARS.hlat = hloc.lat; MARS.hlon = hloc.lon;
E = referenceEllipsoid('wgs84');
dbars = [33.8 -125.5]; dbardist = 100000;
[latout,lonout] = reckon(dbars(1),dbars(2),dbardist,90,E);
    
close all;
% Panel a
figure(1); clf; set(gcf,'position',[200 200 600 300],'color','w');
colormap(flipud(brewermap(128,'Spectral')))
XL = [.5 12.5];  cax = [60 130]; 
monamz = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ordr = [5:12 1:4]; monamz = monamz(ordr);
dg = [.5 .5 .5]; gry = [.7 .7 .7]; ms = 5; fs = 11; XL = [.5 12.5];
set(groot,'DefaultAxesFontSize',11);
P0 = [-.02 0 .5 .9]; P1 = [.5 .4 .48 .55]; P2 = [.5 .1 .48 .285];
cbp = [.08 .84 .3 .02];
AX = [-126.5 -119.5 33.2 39.9];
axes('position',P0);
axesm('MapProjection','Mercator','MapLatLimit',AX([3 4]),'MapLonLimit',AX([1 2]),...
    'PlineLocation',2,'MlineLocation',2,'MeridianLabel','on','ParallelLabel','on',...
    'fontsize',fs,'MlabelParallel','south','GColor','k','FEdgeColor',[.5 .5 .5],'FontColor','k','FLineWidth',1);
framem; gridm; axis off
geoshow(RL.glat,RL.glon,RL.grrl,'DisplayType','texturemap'); caxis(cax); 
geoshow(thiscoast,'FaceColor',gry,'EdgeColor','k');
plotm(dbars(1)+[0 0]-.02,[dbars(2) lonout],'k','linewidth',4);
textm(dbars(1)-.25,mean([dbars(2) lonout]),[int2str(dbardist/1000) ' km'],'horizontalalignment','center','fontsize',fs);
textm(36.95,-121.79,'Monterey','fontsize',fs,'fontweight','bold'); 
textm(36.68,-121.79,'Bay','fontsize',fs,'fontweight','bold');
plotm(MARS.hlat,MARS.hlon,'ko','markersize',5);
cb = linspace(cax(1),cax(2),100); 
axes('position',cbp); 
imagesc(cb,[0 1],[cb;cb]);
set(gca,'Ytick',[],'Xaxislocation','top','fontsize',fs,'Xlim',[cax(1)+1 cax(2)]); 
xlabel('Received level (dB re 1 \muPa)');
 
%Panel b
axes('position',P1);
plot(CI.n.blue(ordr),'k-o','markersize',ms,'markerfacecolor',[.4 .4 .4],'color',[.4 .4 .4],'linewidth',1); hold on;
plot(CI.dd.blue(ordr),'k-o','markersize',ms,'markerfacecolor',[.8 .8 .8],'color',[.8 .8 .8],'linewidth',1);
plot(CI.d.blue(ordr),'k--o','markersize',ms,'markerfacecolor','w','linewidth',1);
legend('Night','Dusk/Dawn','Day');
set(gca,'Ylim',[1 1.14],'Xtick',[1:12],'Xlim',XL,'Xticklabel',[],'Tickdir','out','box','off','fontsize',fs,'color','none');
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl.*[2 2.5]);
ylabel('Blue whale B call index');
frameax

%Panel c
co = get(gca,'colororder'); mb = co(1,:);
axes('position',P2);
d = [CI.n.blue(:)-1 CI.d.blue(:)-1]; dp.blue = 100*(d(:,1)-d(:,2))./d(:,2);
dp.blue(2:6) = 0;
bar(dp.blue(ordr),'FaceColor',mb,'EdgeColor','none'); hold on;
set(gca,'Ylim',[-15 39]);
set(gca,'Xlim',XL,'Xtick',[1:12],'Xticklabel',monamz,'Tickdir','out','box','off','fontsize',fs);
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2);
ylabel('\DeltaCI (%)');

% Annotation
lx = 9.7; ly = [1 10]; 
dg = [0 0 0];  % set to dark gray or black
plot(lx+[0 0],ly,'k','linewidth',2,'color',dg); 
plot(lx,max(ly),'k^','markersize',7,'color',dg,'markerfacecolor',dg);
text(lx+.25,mean(ly)+6.4,'More song','color',dg,'fontsize',11);
text(lx+.25,mean(ly)-.62,'during night','color',dg,'fontsize',11);
plot(lx+[0 0],-ly,'k','linewidth',2,'color',dg); 
plot(lx,-max(ly),'kv','markersize',7,'color',dg,'markerfacecolor',dg);
text(lx+.25,-mean(ly)+1.5,'More song','color',dg,'fontsize',11);
text(lx+.25,-mean(ly)-4,'during day','color',dg,'fontsize',11);
frameax;

% axes positioning and panel labels
axes('position',[0 0 1 1]); set(gca,'visible','off');
text(.02,.94,'A','fontsize',12,'fontweight','bold');
text(.425,.93,'B','fontsize',12,'fontweight','bold');
text(.425,.36,'C','fontsize',12,'fontweight','bold');