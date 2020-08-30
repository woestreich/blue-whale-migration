% investigating diel mvmt patterns of CATS tags in order to address
% reviewer comment
clear all

catsdir = dir('tag_data/CATS/');
locs=[];
j=1;
for i=4:13
    load(['tag_data/CATS/',catsdir(i).name,'/',catsdir(i).name,'_taginfo.mat'])
    ind = find(~isnan(GPS(:,1)) & tagon);
    k=length(ind)-1;
    locs(j:j+k,1) = DN(ind);
    locs(j:j+k,2) = GPS(ind,1);
    locs(j:j+k,3) = GPS(ind,2);
    j=j+k+1;
end

load 'acoustic_data/MARS_hydrophone_location.mat'
x = locs(:,2) - hloc.lat;
y = locs(:,3) - hloc.lon;
dist = sqrt(x.^2 + y.^2);

[~,el] = SolarAzEl(locs(:,1)-(7/24),hloc.lat,hloc.lon,0);
dist_day = dist(el > 0);
dist_night = dist(el < -12);
