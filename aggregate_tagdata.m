%% Aggregate CATS and TDR10 tag data for call and feeding lunge analyses
% Last update: May 18, 2020

% Aggregates across all deployments within each year to avoid the issue of
% very small windows in certain temporal bins (i.e. night) when there are
% CATS deployments that come off early in the night, or before night.

clear all 

%% 2017
CATSdepl = dir(['tag_data/CATS/','*17*']);

% arrays for calls, lunges, and total hours
taghrs = [0,0,0];
calls = [0,0,0];
lunges = [0,0,0];

% cycle through list of all CATS deployments in 2017 w/ calls (A and/or B) and
% lunges 
for f=1:length(CATSdepl)
    disp([num2str(f),'/',num2str(length(CATSdepl))])

    load(['tag_data/CATS/',CATSdepl(f).name,'/',CATSdepl(f).name,'_taginfo.mat'],'DN','tagon','GPS');
    load(['tag_data/CATS/',CATSdepl(f).name,'/',CATSdepl(f).name,'lunges.mat'],'LungeI');
    
    %GPS info
    gpsind = find(~isnan(GPS(:,1)),1);
    if sum(gpsind > 0)
        lat = GPS(gpsind,1); lon = GPS(gpsind,2);
    else
        lat = 36.7125; lon = -122.1868;
    end
    
    %time info
    t1 = DN(1);
    
    %lunge times with highest confidence (LungeC = 3); all 2017 deployment
    %lunge files contain only confidence level 3.
    lungeDN = DN(LungeI);
    
    %load in calls and call times, filtering for A and B calls only
    calltable = readtable(['tag_data/CATS/',CATSdepl(f).name,'/call_detections.txt'],'Delimiter','\t');
    idx = strcmp(calltable.Unit,'A') | strcmp(calltable.Unit,'B');
    calltable = calltable(idx,:);
    calltimes = t1 + (calltable.BeginTime_s_./(60*60*24));
    
    %gets tagon DN's into 1-min resolution array for later calculation of
    %hours at all solar elevations
    ton = DN(tagon);
    tagx = ton(1):(1/24/60):ton(end);
    
    %solar elevation for call, lunge, and tag on times
    [~,callsol] = SolarAzEl(calltimes+(7/24),lat,lon,0);
    [~,lungesol] = SolarAzEl(lungeDN+(7/24),lat,lon,0);
    [~,tagsol] = SolarAzEl(tagx+(7/24),lat,lon,0);
    
    %tag on hours in each solar elevation bin
    taghrs(1) = taghrs(1) + sum(tagsol>0)/60;
    taghrs(2) = taghrs(2) + sum(tagsol<=-12)/60;
    taghrs(3) = taghrs(3) + sum(tagsol<=0 & tagsol>-12)/60;
    
    %calls in each solar elevation bin
    calls(1) = calls(1) + sum(callsol>0);
    calls(2) = calls(2) + sum(callsol<=-12);
    calls(3) = calls(3) + sum(callsol<=0 & callsol>-12);
    
    %lunges in each solar elevation bin
    lunges(1) = lunges(1) + sum(lungesol>0);
    lunges(2) = lunges(2) + sum(lungesol<=-12);
    lunges(3) = lunges(3) + sum(lungesol<=0 & lungesol>-12);
end

%rates
callrate = calls./taghrs;
lungerate = lunges./taghrs;    

%totals
totalhrs = nansum(taghrs);
totalcalls = nansum(calls);
totallunges = nansum(lunges);
depls = 4;

%save
save('tag_data/CATS_TDR_diel_2017.mat','callrate','lungerate','totalhrs','totalcalls','totallunges','depls')

%% 2018
CATSdepl = dir(['tag_data/CATS/','*18*']);
taghrs = [0,0,0];
calls = [0,0,0];
lunges = [0,0,0];

% cycle through list of all CATS deployments in 2018 w/ calls (A and/or B) and
% lunges 
for f=1:length(CATSdepl)
    disp([num2str(f),'/',num2str(length(CATSdepl))])

    load(['tag_data/CATS/',CATSdepl(f).name,'/',CATSdepl(f).name,' 10Hzprh.mat'],'DN','tagon','GPS');
    load(['tag_data/CATS/',CATSdepl(f).name,'/',CATSdepl(f).name,'lunges.mat'],'LungeI','LungeC');
    
    %GPS info
    gpsind = find(~isnan(GPS(:,1)),1);
    if sum(gpsind > 0)
        lat = GPS(gpsind,1); lon = GPS(gpsind,2);
    else
        lat = 36.7125; lon = -122.1868;
    end
    
    %time info
    t1 = DN(1);
    
    %lunge times with highest confidence (LungeC = 3)
    lungeDN = DN(LungeI);
    lungeDN = lungeDN(LungeC == 3);
    
    %load in calls and call times, filtering for A and B calls only
    calltable = readtable(['tag_data/CATS/',CATSdepl(f).name,'/call_detections.txt'],'Delimiter','\t');
    idx = strcmp(calltable.Unit,'A') | strcmp(calltable.Unit,'B');
    calltable = calltable(idx,:);
    calltimes = t1 + (calltable.BeginTime_s_./(60*60*24));
    
    %gets tagon DN's into 1-min resolution array for later calculation of
    %hours at all solar elevations
    ton = DN(tagon);
    tagx = ton(1):(1/24/60):ton(end);
    
    %solar elevation for call, lunge, and tag on times
    [~,callsol] = SolarAzEl(calltimes+(7/24),lat,lon,0);
    [~,lungesol] = SolarAzEl(lungeDN+(7/24),lat,lon,0);
    [~,tagsol] = SolarAzEl(tagx+(7/24),lat,lon,0);
    
    %tag on hours in each solar elevation bin
    taghrs(1) = taghrs(1) + sum(tagsol>0)/60;
    taghrs(2) = taghrs(2) + sum(tagsol<=-12)/60;
    taghrs(3) = taghrs(3) + sum(tagsol<=0 & tagsol>-12)/60;
    
    %calls in each solar elevation bin
    calls(1) = calls(1) + sum(callsol>0);
    calls(2) = calls(2) + sum(callsol<=-12);
    calls(3) = calls(3) + sum(callsol<=0 & callsol>-12);
    
    %lunges in each solar elevation bin
    lunges(1) = lunges(1) + sum(lungesol>0);
    lunges(2) = lunges(2) + sum(lungesol<=-12);
    lunges(3) = lunges(3) + sum(lungesol<=0 & lungesol>-12);
end

%load in 2018 TDR w/ calls
disp('TDR')
load('tag_data/TDR10/Bm181021-TDR11/Bm181021-TDR11_taginfo.mat','DN','tagon');
load('tag_data/TDR10/Bm181021-TDR11/Bm181021-TDR11 lunges.mat','LungeDN','LungeC');
load('tag_data/TDR10/Bm181021-TDR11/Bm181021-TDR11_32HzAcc_calls.mat','callDN');
gps = xlsread('tag_data/TDR10/Bm181021-TDR11/Bm181021-TDR11 GPS.xlsx');

%subset TDR variables of interest to foraging period
callDN = callDN(callDN < datenum('1 Nov 2018'));
idx = DN < datenum('1 Nov 2018');
DN = DN(idx); 
tagon = tagon(idx);

%subset lunges to only highest confidence level (3)
LungeDN = LungeDN(LungeC == 3);

%gets tagon DN's into 1-min resolution array for later calculation of
%hours at all solar elevations
ton = DN(tagon);
tagx = ton(1):(1/24/60):ton(end);

%get nearest lat and lon for all calls, lunges, and tagon times (for solar
%elevation calculations
clat = []; clon = [];
llat = []; llon = [];
tlat = []; tlon = [];
for j=1:numel(callDN)
    cdiff = callDN(j) - gps(:,6);
    [~,idx] = min(abs(cdiff));
    clat(j) = gps(idx,2);
    clon(j) = gps(idx,3);
end
for j=1:numel(LungeDN)
    ldiff = LungeDN(j) - gps(:,6);
    [~,idx] = min(abs(ldiff));
    llat(j) = gps(idx,2);
    llon(j) = gps(idx,3);
end
for j=1:numel(tagx)
    tdiff = tagx(j) - gps(:,6);
    [~,idx] = min(abs(tdiff));
    tlat(j) = gps(idx,2);
    tlon(j) = gps(idx,3);
end

%solar elevation for call, lunge, and tag on times
[~,callsol] = SolarAzEl(callDN+(7/24),clat,clon,0);
[~,lungesol] = SolarAzEl(LungeDN+(7/24),llat,llon,0);
[~,tagsol] = SolarAzEl(tagx+(7/24),tlat,tlon,0);

%tag on hours in each solar elevation bin
taghrs(1) = taghrs(1) + sum(tagsol>0)/60;
taghrs(2) = taghrs(2) + sum(tagsol<=-12)/60;
taghrs(3) = taghrs(3) + sum(tagsol<=0 & tagsol>-12)/60;
    
%calls in each solar elevation bin
calls(1) = calls(1) + sum(callsol>0);
calls(2) = calls(2) + sum(callsol<=-12);
calls(3) = calls(3) + sum(callsol<=0 & callsol>-12);
    
%lunges in each solar elevation bin
lunges(1) = lunges(1) + sum(lungesol>0);
lunges(2) = lunges(2) + sum(lungesol<=-12);
lunges(3) = lunges(3) + sum(lungesol<=0 & lungesol>-12);

%NOW FINAL RATE AND TOTAL CALCULATIONS W/ CATS + TDR
%rates
callrate = calls./taghrs;
lungerate = lunges./taghrs;    

%totals
totalhrs = nansum(taghrs);
totalcalls = nansum(calls);
totallunges = nansum(lunges);
depls = 6;

%save
save('CATS_TDR_diel_2018.mat','callrate','lungerate','totalhrs','totalcalls','totallunges','depls')

%% 2019
CATSdepl = dir(['tag_data/CATS/','*19*']);
taghrs = [0,0,0];
calls = [0,0,0];
lunges = [0,0,0];

% cycle through list of all CATS deployments in 2019 w/ calls (A and/or B) and
% lunges 
for f=1:length(CATSdepl)
    disp([num2str(f),'/',num2str(length(CATSdepl))])
    
    load(['tag_data/CATS/',CATSdepl(f).name,'/',CATSdepl(f).name,' 10Hzprh.mat'],'DN','tagon','GPS');
    load(['tag_data/CATS/',CATSdepl(f).name,'/',CATSdepl(f).name,'lunges.mat'],'LungeI','LungeC');
    
    %GPS info
    gpsind = find(~isnan(GPS(:,1)),1);
    if sum(gpsind > 0)
        lat = GPS(gpsind,1); lon = GPS(gpsind,2);
    else
        lat = 36.7125; lon = -122.1868;
    end
    
    %time info
    t1 = DN(1);
    
    %lunge times with highest confidence (LungeC = 3)
    lungeDN = DN(LungeI);
    lungeDN = lungeDN(LungeC == 3);
    
    %load in calls and call times, filtering for A and B calls only
    calltable = readtable(['tag_data/CATS/',CATSdepl(f).name,'/call_detections.txt'],'Delimiter','\t');
    idx = strcmp(calltable.Unit,'A') | strcmp(calltable.Unit,'B');
    calltable = calltable(idx,:);
    calltimes = t1 + (calltable.BeginTime_s_./(60*60*24));
    
    %gets tagon DN's into 1-min resolution array for later calculation of
    %hours at all solar elevations
    ton = DN(tagon);
    tagx = ton(1):(1/24/60):ton(end);
    
    %solar elevation for call, lunge, and tag on times
    [~,callsol] = SolarAzEl(calltimes+(7/24),lat,lon,0);
    [~,lungesol] = SolarAzEl(lungeDN+(7/24),lat,lon,0);
    [~,tagsol] = SolarAzEl(tagx+(7/24),lat,lon,0);
    
    %tag on hours in each solar elevation bin
    taghrs(1) = taghrs(1) + sum(tagsol>0)/60;
    taghrs(2) = taghrs(2) + sum(tagsol<=-12)/60;
    taghrs(3) = taghrs(3) + sum(tagsol<=0 & tagsol>-12)/60;
    
    %calls in each solar elevation bin
    calls(1) = calls(1) + sum(callsol>0);
    calls(2) = calls(2) + sum(callsol<=-12);
    calls(3) = calls(3) + sum(callsol<=0 & callsol>-12);
    
    %lunges in each solar elevation bin
    lunges(1) = lunges(1) + sum(lungesol>0);
    lunges(2) = lunges(2) + sum(lungesol<=-12);
    lunges(3) = lunges(3) + sum(lungesol<=0 & lungesol>-12);
end

%load in 2019 TDR w/ calls
disp('TDR')
load('tag_data/TDR10/Bm190916-TDR14/Bm190916-TDR14_taginfo.mat','DN','tagon');
load('tag_data/TDR10/Bm190916-TDR14/Bm190916-TDR14 lunges.mat','LungeDN','LungeC');
load('tag_data/TDR10/Bm190916-TDR14/Bm190916-TDR14_32HzAcc_calls.mat','callDN');
gps = xlsread('tag_data/TDR10/Bm190916-TDR14/Bm190916-TDR14 GPS.xlsx');

%subset TDR variables of interest to foraging period
callDN = callDN(callDN < datenum('23 Sept 2019'));
idx = DN < datenum('23 Sept 2019');
DN = DN(idx); 
tagon = tagon(idx);

%gets tagon DN's into 1-min resolution array for later calculation of
%hours at all solar elevations
ton = DN(tagon);
tagx = ton(1):(1/24/60):ton(end);

%subset lunges to only highest confidence level (3)
LungeDN = LungeDN(LungeC == 3);

%get nearest lat and lon for all calls, lunges, and tagon times (for solar
%elevation calculations
clat = []; clon = [];
llat = []; llon = [];
tlat = []; tlon = [];
for j=1:numel(callDN)
    cdiff = callDN(j) - gps(:,6);
    [~,idx] = min(abs(cdiff));
    clat(j) = gps(idx,2);
    clon(j) = gps(idx,3);
end
for j=1:numel(LungeDN)
    ldiff = LungeDN(j) - gps(:,6);
    [~,idx] = min(abs(ldiff));
    llat(j) = gps(idx,2);
    llon(j) = gps(idx,3);
end
for j=1:numel(tagx)
    tdiff = tagx(j) - gps(:,6);
    [~,idx] = min(abs(tdiff));
    tlat(j) = gps(idx,2);
    tlon(j) = gps(idx,3);
end

%solar elevation for call, lunge, and tag on times
[~,callsol] = SolarAzEl(callDN+(7/24),clat,clon,0);
[~,lungesol] = SolarAzEl(LungeDN+(7/24),llat,llon,0);
[~,tagsol] = SolarAzEl(tagx+(7/24),tlat,tlon,0);

%tag on hours in each solar elevation bin
taghrs(1) = taghrs(1) + sum(tagsol>0)/60;
taghrs(2) = taghrs(2) + sum(tagsol<=-12)/60;
taghrs(3) = taghrs(3) + sum(tagsol<=0 & tagsol>-12)/60;
    
%calls in each solar elevation bin
calls(1) = calls(1) + sum(callsol>0);
calls(2) = calls(2) + sum(callsol<=-12);
calls(3) = calls(3) + sum(callsol<=0 & callsol>-12);
    
%lunges in each solar elevation bin
lunges(1) = lunges(1) + sum(lungesol>0);
lunges(2) = lunges(2) + sum(lungesol<=-12);
lunges(3) = lunges(3) + sum(lungesol<=0 & lungesol>-12);

%NOW FINAL RATE AND TOTAL CALCULATIONS W/ CATS + TDR
%rates
callrate = calls./taghrs;
lungerate = lunges./taghrs;    

%totals
totalhrs = nansum(taghrs);
totalcalls = nansum(calls);
totallunges = nansum(lunges);
depls = 3;

%save
save('CATS_TDR_diel_2019.mat','callrate','lungerate','totalhrs','totalcalls','totallunges','depls')
