%% Generate Figure 1 
% Panel a: example blue whale song 
% Panel b: spectral statistics, CI calculation frequenices, and reference
% map
% Last update: May 18, 2020

clear all; close all;

% example blue whale song
load 'acoustic_data/bluewhalesong_ltsa_20161101h1.mat'; 
load 'acoustic_data/bluewhalesong_xl.mat'; 

% Aug-Dec (peak blue whale song season) 50th and 90th percentils
load 'acoustic_data/AugDecPercentiles.mat';

% land areas for map
landareas = shaperead('landareas.shp','UseGeoCoords',true);
load 'acoustic_data/MARS_hydrophone_location.mat'
gry = .7+[0 0 0]; dg = .5+[0 0 0];

% set up figure and positioning
figure(1); clf; 
set(gcf,'position',[200 200 650 350],'color','w');
fs = 11;
set(groot,'DefaultAxesFontSize',fs);
P1 = [.06 .125 .67 .72];
cbp = [P1(1)+.04 P1(2)+P1(4)+.02 P1(3)-.08 .015];
P2 = [P1(1)+P1(3)+.01 P1(2) .24 P1(4)];
insetP = [.745 .535 .3 .3]; gry = [.7 .7 .7];
cmp = cmocean('ice',256); 
colormap(cmp); 
fr = [6 93]; yl = fr+[1 -1];
cax = [73 105]; clevs = [0 linspace(cax(1),cax(2),15) 200];

% LTSA
axes('position',P1);
fidx = find(L.freq >= fr(1) & L.freq < fr(2));
tidx = find(L.dn >= xl(1) & L.dn <= xl(2));
wtime = [L.time(tidx) - min(L.time(tidx))]/60; 
[cc,h] = contourf(wtime,(L.freq(fidx)),L.ltsa(fidx,tidx),clevs); set(h(:),'edgecolor','none'); 
axis tight; caxis(cax); ylabel('Frequency (Hz)'); xlabel('Time (minutes)');
set(gca,'Ylim',yl,'Xlim',[0 max(wtime)+1.25],'color',[.7 .7 .7],'fontsize',fs,'box','off');
uf = [10.5 14.5+[0:4]*14.5 81]; un = {'C','B','B_2','B_3','B_4','B_5','A'};
for u = 1:numel(uf);
    text(max(wtime)+.15,uf(u),un{u},'fontsize',11); 
end
frameax;
cb = linspace(cax(1),cax(2),100);
axes('position',cbp);
imagesc(cb,[0 1],[cb;cb]); 
set(gca,'Ytick',[],'Xaxislocation','top','fontsize',fs);
xlabel('Spectrum Level (dB re 1 \muPa^2 Hz^{-1})'); 

% Spectrum level statistics
dg = [.5 .5 .5];
axes('position',P2);
fidx = find(Q.freq >= fr(1) & Q.freq < fr(2));
plot(Q.pctile50(fidx),Q.freq(fidx),'k','linewidth',1); 
axis tight; 
set(gca,'Xlim',[70 100])
hold on;

% algorithm frequencies
cif.blue = [37 43 44 50]; axl = get(gca,'Xlim');
for F = 1:numel(cif.blue);
    [x,ia,ib] = intersect(cif.blue(F),Q.freq);
    plot([Q.pctile50(ib)+.2 axl(2)],cif.blue(F)+[0 0],'k','color',dg);
end

set(gca,'Ylim',yl,'box','off','Ytick',[],'fontsize',fs,'Tickdir','out');
frameax
axes('position',P2);
plot(Q.pctile90(fidx),Q.freq(fidx),'k--','linewidth',1); 
axis tight;
set(gca,'Xlim',[70 100])
xlabel('dB re 1 \muPa^2 Hz^{-1}'); 
set(gca,'Ylim',yl,'box','off','Yaxislocation','right','color','none','Ytick',[],'fontsize',fs,'Tickdir','out');

% reference map
addpath('m_map/')
axes('position',insetP);
usemmap = 1;
if usemmap
    m_proj('ortho','lat',42','long',-150');    
    m_plot(hloc.lon,hloc.lat,'ko','markerfacecolor','g','markersize',5,'color','none');
    m_coast('patch',gry,'Edgecolor','none'); hold on;
    m_grid('linest','-','xticklabels',[],'yticklabels',[]);
    set(gca,'color','none');
else
    axesm('ortho','Origin',[42 -150 0],'Grid','on','Glinestyle','-','Gcolor',dg,'Frame','off'); %,'MlineLocation',60,'PlineLocation',30);
    geoshow(landareas,'FaceColor',gry,'EdgeColor','none');
    plotm(hloc.lat,hloc.lon,'ko','markerfacecolor','g','markersize',5);
end

axes('position',[0 0 1 1]); set(gca,'visible','off');
lx = P2(1)+[.01 .05]; ly = .2;
plot(lx-.005,ly+[0 0],'k--','linewidth',1); text(lx(2)+.005,ly+0.015,'90th','fontsize',fs,'fontweight','normal'); text(lx(2)+.005,ly-0.015,'pctl','fontsize',fs,'fontweight','normal');
hold on;
ly = ly+.08;
plot(lx-.005,ly+[0 0],'k','linewidth',1); text(lx(2)+.005,ly+0.015,'50th','fontsize',fs,'fontweight','normal'); text(lx(2)+.005,ly-0.015,'pctl','fontsize',fs,'fontweight','normal');
text(P2(1)+.5*P2(3),P2(2)+P2(4)+.04,['August ' char(8212) ' December'],'horizontalalignment','center','fontsize',12);
axis([0 1 0 1]); set(gca,'visible','off');
text(P1(1),.97,'a','fontweight','bold','fontsize',12);
text(P2(1),.97,'b','fontweight','bold','fontsize',12);
text(.87,.705,'North','fontsize',11);
text(.87,.675,'Pacific','fontsize',11);
