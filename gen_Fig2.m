%% Generate Figure 2
% Panel A: modeled received level results for MARS hydrophone
%
% Panel B: average annual cycle of call index (CI)
%
% Panel C: average annual cycle of CI night:day ratio
%
% First section (Lines X-Y) calculates CI, CI ratio, and quartiles + 
% percentiles for presentation in panels B and C. This section also
% generates 4 exploratory figures for CI and CI ratio.
%
% Second section (Lines X-Y) saves a csv to pass to R (where we conduct 
% 2-sided t-tests to test for month-to-month significant changes in CI and
% CI ratio (also presented on panels B and C).
%
% Third (final) section (Lines X-Y) creates 3-panel figure
%
% Last update: August 30, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; 

% Daily power spectral density binned by solar elevation 
load 'acoustic_data/MARS_SpectrumLevel_Daily_bySolarElevation.mat'
    
% Handle for month
dv = datevec(D.time); mo = dv(:,2); 

% Quartiles for daily CI and daily categorical CI
[nfreq,ntime,ns] = size(D.sm);
% First examine counts, which should be independent of frequency
sn = {'n','dd','d','all'};
pcount = squeeze(D.ct(1,:,:)); pcount = pcount';
pcount(4,:) = sum(pcount(1:3,:));
        
figure(1); clf; 

mincount = [380 100 500 1200]; % night, dd, day, total
        
for S = 1:4;
    eval(['subplot(41' int2str(S) ')']);
    plot(D.time,pcount(S,:),'-o');
    ylabel('Observation count');
    axis tight; datetick('x'); hold on;
    xl = get(gca,'Xlim'); plot(xl,mincount(S)+[0 0],'r--');
    idx = find(pcount(S,:) < mincount(S));
    spct = round(10000*numel(idx)/ntime)/100;
    title([sn{S} ': ' num2str(spct) '% of days < minimum threshold']);
end

% Get daily call indices by SE category, screened
clear DCI; DCI.time = D.time; DCI.sn = sn; 
for S = 1:3;
    sm = D.sm(:,:,S); ct = D.ct(:,:,S);
    % remove data for days with insufficient sampling
    xcl = find(ct(1,:) < mincount(S)); ct(:,xcl) = NaN;
    % Place into structure
    clear L; L.time = D.time; L.freq = D.freq; L.ltsa = sm./ct;
    % Get the daily call index for this SE category
    c = call_index(L); DCI.blue(S,:) = c.blue;
end
% Get daily call index for all data, regardless of SE category
sm = sum(D.sm,3); ct = sum(D.ct,3);
S = 4;
xcl = find(ct(1,:) < mincount(S)); ct(:,xcl) = NaN;
clear L; L.time = D.time; L.freq = D.freq; L.ltsa = sm./ct;
c = call_index(L); DCI.blue(S,:) = c.blue;
q = DCI.blue([1 3],:); % isolate night and day
q = q - min(q(:));  % scale range above minimum
S = 5; DCI.blue(S,:) = q(1,:)./q(2,:);

% Plot the daily values by category
sn = {'n','dd','d','all','ratio'};

figure(2); clf; 
for S = 1:5;
    eval(['subplot(23' int2str(S) ')']);
    plot(DCI.time,DCI.blue(S,:),'o');
    ng = numel(find(~isnan(DCI.blue(S,:))));
    title([sn{S} ': ' int2str(ng) ' good indices']);
    datetick('x'); axis tight;
end

% To store quartiles by month for each of 3 SE categories, plus
% overall CI and CI ratio, use a 3D matrix.  The result from a
% single month/SE analysis is a row vector of size [1 x 3].  
% So, Q has dimensions {month x quartile x SE bin / overall];
% Also store Q_mean (monthly mean values for CI and delta CI)
        
dv = datevec(DCI.time);
clear Q  % Quartiles for daily (all) CI and delta CI; SE categories
pctl = [.1 .25 .5 .75 .9]; % percentiles
for M = 1:12;
    k = find(dv(:,2) == M);
    % Quartiles by SE category and overall
    for S = 1:5;
        c = DCI.blue(S,k);
        c(isnan(c)) = []; c = sort(c);
        qi = round(pctl*length(c)); % quartile indices
        Q(M,:,S) = c(qi); % quartiles (+ 10th & 90th pctls)
    end
    for S = 4:5;
        c = DCI.blue(S,k); c(isnan(c)) = []; 
        Q_mean(M,S-3)=mean(c); 
    end
end
        
% Quick look
figure(3); clf; 
ttls = {'Night','Dusk/Dawn','Day','All','Delta CI'};
bufr = .3;
for S = 1:5;
    eval(['subplot(23' int2str(S) ')']);
    c = Q(:,:,S); % c = c(:);
    for M = 1:12;
        cq = c(M,:);
        by = cq([2 2 4 4 2]); bx = M+[-1 1 1 -1 -1]*bufr;
        plot(bx,by,'k'); hold on;
        plot(M+[-1 1]*bufr, cq(3)+[0 0],'k');
        if S < 5; set(gca,'Ylim',[.98 1.2]); end
    end
end

figure(4); close;
ordr = [5:12 1:4]; % Month order
fs = 11; set(groot,'DefaultAxesFontSize',fs);
lw = 1; lcol = [.5 .5 .5]; bcol = [.3 .3 .3];

figure(4); clf;
set(gcf,'position',[200 200 300 300]);
px = .175; pw = .79; ph = .43; py = .98-[1 2]*ph - [0 .015];
P1 = [px py(1) pw ph]; P2 = [px py(2) pw ph];
dg = [.5 .5 .5]; gry = [.7 .7 .7]; ms = 5; 
XL = [.5 12.5]; bufr = .3; 

monamz = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
monamz = monamz(ordr);
XT = [.5:1:11.5];

% CI
axes('position',P1);
c = Q(:,:,4); % overall CI
c = c(ordr,:); % ordered May-Dec_Jan-Apr
for M = 1:12;
    cq = c(M,:);
    by = cq([2 2 4 4 2]); bx = M+[-1 1 1 -1 -1]*bufr;
    plot(bx,by,'k','linewidth',lw,'color',lcol); hold on;
    plot(M+[-1 1]*bufr, cq(3)+[0 0],'k','linewidth',lw,'color',lcol);
end
set(gca,'Ylim',[.99 1.165]);
hold on; plot(XL,[1 1],'k--');
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2); 
set(gca,'Xlim',XL,'Xtick',XT,'Xticklabel',[],'Xgrid','on',...
    'Tickdir','out','box','off','fontsize',fs);
ylabel('CI'); 
frameax;
 
% CI ratio
axes('position',P2);
c = Q(:,:,5); % CI ratio
c(2:6,:) = NaN;
c = c(ordr,:); % ordered May-Dec_Jan-Apr
for M = 1:12;
    cq = c(M,:);
    if ~all(isnan(cq))
        by = cq([1 1 3 3 1]); bx = M+[-1 1 1 -1 -1]*bufr;
        plot(bx,by,'k','linewidth',lw,'color',lcol); hold on;
        plot(M+[-1 1]*bufr, cq(2)+[0 0],'k','linewidth',lw,'color',lcol);
    end
end
set(gca,'Ylim',[0 2.25]);
hold on; plot(XL,[1 1],'k--');
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2); 
set(gca,'Xlim',XL,'Xtick',XT,'Xticklabel',[],'Xgrid','on',...
    'Tickdir','out','box','off','fontsize',fs);
ylabel('CI_{night} : CI_{day}'); 
for M = 1:12;
    text(XT(M)+.5,-.24,monamz{M},'Horizontalalignment','center','rotation',90,'fontsize',fs);
end
 
frameax;
%% csv for statistics in R
% save csv for use in R
ci_daily = [DCI.blue(4,:); DCI.blue(5,:); dv(:,2)']';
csvwrite('ci_daily.csv',ci_daily);    
        
%% Generate Figure 2 for main text     
% Received level model output
RL = load('acoustic_data/RL_Sc01_November.mat');
% Adjustment to RL (TL was modeled with SL=186 dB; Thode et al. suggest 
% SL=171 dB). Because RL = SL - TL, RL can be adjusted as follows: 
RL.grrl = RL.grrl - 15;

% West coast for map
thiscoast = shaperead('thiscoast', 'UseGeoCoords', true);

% MARS hydrophone location
load acoustic_data/MARS_hydrophone_location; MARS.hlat = hloc.lat; MARS.hlon = hloc.lon;
E = referenceEllipsoid('wgs84');
%dbars = [33.8 -125.5]; 
dbars = [35.3 -123.9]; 
dbardist = 100000;
[latout,lonout] = reckon(dbars(1),dbars(2),dbardist,90,E);

close all;

% Panel a
figure(1); clf; set(gcf,'position',[200 200 600 300],'color','w');
cmap=flipud(brewermap(128,'Spectral'));
cmap(1,:) = [1 1 1];
colormap(cmap)
XL = [.5 12.5];  cax = [78 120]; 
monamz = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ordr = [5:12 1:4]; monamz = monamz(ordr);
dg = [.5 .5 .5]; gry = [.7 .7 .7]; ms = 5; fs = 11; XL = [.5 12.5];
set(groot,'DefaultAxesFontSize',11);
P0 = [.03 .1 .4 .7]; P1 = [.5 .5 .48 .35]; P2 = [.5 .14 .48 .35];
cbp = [.08 .84 .3 .02];
AX = [-124.7 -120.8 35.05 37.95];
axes('position',P0);
axesm('MapProjection','Mercator','MapLatLimit',AX([3 4]),'MapLonLimit',AX([1 2]),...
    'PlineLocation',1,'MlineLocation',2,'MeridianLabel','on','ParallelLabel','on',...
    'fontsize',fs,'MlabelParallel','south','GColor','k','FEdgeColor',[.5 .5 .5],'FontColor','k','FLineWidth',1);
framem; gridm; axis off
geoshow(RL.glat,RL.glon,RL.grrl,'DisplayType','texturemap'); caxis(cax); 
geoshow(thiscoast,'FaceColor',gry,'EdgeColor','k');
plotm(dbars(1)+[0 0]-.02,[dbars(2) lonout],'k','linewidth',4);
textm(dbars(1)-.15,mean([dbars(2) lonout]),[int2str(dbardist/1000) ' km'],'horizontalalignment','center','fontsize',fs);
textm(36.85,-121.79,'Monterey','fontsize',fs,'fontweight','bold'); 
textm(36.63,-121.79,'Bay','fontsize',fs,'fontweight','bold');
plotm(MARS.hlat,MARS.hlon,'ko','markersize',5);

axes('position',P0);
axesm('MapProjection','Mercator','MapLatLimit',AX([3 4]),'MapLonLimit',AX([1 2]),...
    'MeridianLabel','off','ParallelLabel','off','FEdgeColor',[.5 .5 .5],'FLineWidth',1);
framem; axis off

cb = linspace(cax(1),cax(2),100); 
axes('position',cbp); 
imagesc(cb,[0 1],[cb;cb]);
set(gca,'Ytick',[],'Xaxislocation','top','fontsize',fs,'Xlim',[cax(1)+1 cax(2)]); 
xlabel('Received level (dB re 1 \muPa)');
 
%Panel b
ci_pvs = dlmread('ci_pvs.csv',',',1,0);
axes('position',P1);
c = Q(:,:,4); % overall CI
c = c(ordr,:); % ordered May-Dec_Jan-Apr
c_mean = Q_mean(:,1); % mean CI; % ordered May-Dec_Jan-Apr
c_mean = c_mean(ordr,:);
ci_pvs = ci_pvs(ordr,:);
for M = 1:12;
    cq = c(M,:);
    by = cq([2 2 4 4 2]); bx = M+[-1 1 1 -1 -1]*bufr;
    rectangle('Position',[M-.2,cq(1),.4,cq(5)-cq(1)],'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]); hold on;
    plot(bx,by,'k','linewidth',lw,'color','k'); hold on;
    plot(M+[-1 1]*bufr, cq(3)+[0 0],'k','linewidth',lw,'color','k');
    plot(M, c_mean(M), 'ko', 'MarkerFaceColor','k','MarkerSize',4)
    %indicate significant monthly transitions
    if ci_pvs(M,1) < 0.05
        plot(M+0.5,1,'k*')
    end
end
set(gca,'Ylim',[.99 1.21]);
hold on; plot(XL,[1 1],'k--');
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2); 
set(gca,'Xlim',XL,'Xtick',XT,'Xticklabel',[],'Xgrid','on',...
    'Tickdir','out','box','off','fontsize',fs);
ylabel('CI'); 
frameax;

%Panel c
co = get(gca,'colororder'); mb = co(1,:);
axes('position',P2);
c = Q(:,:,5); % Delta CI
c(2:6,:) = NaN;
c = c(ordr,:); % ordered May-Dec_Jan-Apr
c_mean = Q_mean(:,2);
c_mean(2:6) = NaN;
c_mean = c_mean(ordr,:);
for M = 1:12;
    cq = c(M,:);
    if ~all(isnan(cq))
        by = cq([2 2 4 4 2]); bx = M+[-1 1 1 -1 -1]*bufr;
        rectangle('Position',[M-.2,cq(1),.4,cq(5)-cq(1)],'FaceColor',[.8 .8 .8],'EdgeColor',[.8 .8 .8]); hold on;
        plot(bx,by,'k','linewidth',lw,'color','k'); 
        plot(M+[-1 1]*bufr, cq(3)+[0 0],'k','linewidth',lw,'color','k');
        plot(M,c_mean(M),'ko','MarkerFaceColor','k','Markersize',4)
    end
    %indicate significant monthly transitions
    if (ci_pvs(M,2) < 0.05) & (M > 2) & (M < 9) 
        plot(M+0.5,1,'k*')
    end
end
set(gca,'Ylim',[0.6 1.98]);
hold on; plot(XL,[1 1],'k--');
tl = get(gca,'Ticklength'); set(gca,'Ticklength',tl*2); 
set(gca,'Xlim',XL,'Xtick',XT,'Xticklabel',[],'Xgrid','on',...
    'Tickdir','out','box','off','fontsize',fs);
ylabel('CI_{night} : CI_{day}'); 
for M = 1:12;
    text(XT(M)+.5,.4,monamz{M},'Horizontalalignment','center','rotation',90,'fontsize',fs);
end

frameax;

% axes positioning and panel labels
axes('position',[0 0 1 1]); set(gca,'visible','off');
text(.02,.94,'A','fontsize',12,'fontweight','bold');
text(.42,.93,'B','fontsize',12,'fontweight','bold');
text(.42,.39,'C','fontsize',12,'fontweight','bold');