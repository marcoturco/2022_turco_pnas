%% start
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_frap=1950:2021;
years_study=1971:2021;
years_mtbs=1984:2020;
years_westerling=1980:2004;
years_williams=1972:2018;
dir_fire='~/Dropbox/estcena/scripts/fires_california/data_def/fires/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
dir_westerling='~/Documents/virtualBox/lavori/drought_fire_cal/';
dir_shp='~/Dropbox/estcena/scripts/fires_california/data_def/shapefiles/baileys/';
dir_frap = '~/Documents/dati/fire_us/';
font_size=12;

%% FRAP data
% one way to 'map' fire polygons
nomefile = [dir_fire,'forestlayer.mat'];
load(nomefile)

filename = [dir_frap, 'fire21_sierra_northcoast.shp'];
frap = shaperead(filename, 'UseGeoCoords', true); % imports lat/lon data

% create a grid that matches the forest layer
f=find(lat>=34.5 & lat<=42);
f2=find(lon<-118);
lat=lat(f);lon=lon(f2);
forest=forest(f,f2);
numyears=72; % 1950-2021
inityear=1950;
blankgrid=zeros(length(f),length(f2),numyears*12);
[lon,lat]=meshgrid(lon,lat);
% loop over each FRAP fire; I have frap as a structure in MATLAB


% find(frap.ALARM_DATE == '' )
% aux=char(frap(:).ALARM_DATE);
% find(aux == '')


for i=1:length(frap)
    i
    fd=datevec(frap(i).ALARM_DATE);
    
    
    if ~isempty(fd)
        if fd(1)>=inityear & fd<=inityear+numyears
            % can bin by month if needed, here I am just doing all fires in a year
            kk=find(isnan(frap(i).Lon));
            ilast=0;
            for j=1:length(kk)
                p=inpolygon(lon,lat,frap(i).Lon(ilast+1:kk(j)-1),frap(i).Lat(ilast+1:kk(j)-1));
                blankgrid(:,:,(fd(1)-inityear)*12+fd(2))=blankgrid(:,:,(fd(1)-inityear)*12+fd(2))+p;
                ilast=kk(j);
            end
        end
    end
end

%contourm(lat,lon,double(forest))
%
contourm(lat,lon,double(blankgrid(:,:,600-4)))

% you can then mask by region and calculate burned area by scaling by the
% area of the grids

R = earth_radius(lat,lon);
A = cdtarea(lat,lon,'km2');

blankgrid_domain=blankgrid*NaN;
frap_forest_sierra_ncoast=zeros(numyears*12,1)*NaN;
for i=1:size(blankgrid,3)
    frap_forest_sierra_ncoast(i)=sum(sum(double(blankgrid(:,:,i)).*A,1),2);
end

plot(frap_forest_sierra_ncoast)

frap_forest_sierra_ncoast_year=zeros(numyears,1)*NaN;
for iyear=1:length(years_frap) 
  i1 = (iyear - 1) * 12 + 5;
  i2 = (iyear - 1) * 12 + 9;
  frap_forest_sierra_ncoast_year(iyear)=sum(frap_forest_sierra_ncoast(i1:i2));
end


figure; plot(frap_forest_sierra_ncoast_year)
FIRE=frap_forest_sierra_ncoast_year;
savename = [dir_fire,'frap_forest_sierra_ncoast_year_1950_2021.mat'];
save(savename,'FIRE')
[C,IA,IB] = intersect(years_study,years_frap);
FIRE=frap_forest_sierra_ncoast_year(IB);
%FIRE=frap_forest_sierra_ncoast_year;
savename = [dir_fire,'frap_forest_sierra_ncoast_year.mat'];
save(savename,'FIRE')


%% compare other fire datasets


namefile = [dir_fire,'frap_forest_sierra_ncoast_year.mat'];
load(namefile,'FIRE')


% all frap
nomefile = [dir_fire,'all_fire21_sierra_northcoast.mat'];
load(nomefile,'BA')
BA=BA';
for iyear=1:length(years_study)
    i1=(iyear-1)*12+5
    i2=i1+4
    dum(iyear)=nansum(BA(i1:i2));
end
FRAP=dum'*0.0040468564224; %acres to km2;

figure; hold on; 
plot(FIRE)
plot(FRAP,'r')
sum(FIRE)/sum(FRAP)

% mtbs
nomefile = [dir_fire,'mtbs21_sierra_northcoast.mat'];
load(nomefile,'BA')
[C,IA,IB] = intersect(years_study,years_mtbs);
BA=BA';
for iyear=1:length(years_mtbs)
    i1=(iyear-1)*12+5
    i2=i1+4
    dum(iyear)=nansum(BA(i1:i2));
end
MTBS=zeros(length(years_study),1)*NaN;
MTBS(IA)=dum(IB); clear BAS

% westerling
filename = [dir_westerling,'ALL.ACRE.8004.xlsx']; %
[NUM,TXT,RAW]=xlsread(filename);
clear TXT RAW
coord(:,1)=NUM(2,3:end)';
coord(:,2)=NUM(1,3:end)';
BAm_all=NUM(3:end,3:end)*(4.04685642/1000); %acres to kmq
% load california
filename = [dir_shp,'Baileys_ecoregions_cal_wgs84.shp'];
S = shaperead(filename);
%geoshow(S)
for izone=1:size(S,1)
    dum=S(izone).X';
    xv=[dum;dum(1)]; clear dum
    dum=S(izone).Y';
    yv=[dum;dum(1)]; clear dum
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    %figure;drawStations(coord(in,:))
    % monthly aggregation
    k=0;
    for iyear=years_westerling
        for im=1:12
            k=k+1;
            BAm(k,izone)=nansum(BAm_all(k,in),2);
        end
    end
end
% MJJAS aggregation
BAs=zeros(length(years_westerling),size(S,1));
for izone=1:size(S,1)
    for iyear=1:length(years_westerling)
        i1=(iyear-1)*12+5
        i2=i1+4
        BAs(iyear,izone)=nansum(BAm(i1:i2,izone));
    end
end
[C,IA,IB] = intersect(years_study,years_westerling);
WEST=zeros(length(years_study),1)*NaN;
WEST(IA)=BAs(IB); clear BAs

%williams et al. (2019)
nomefile = [dir_fire,'forestburned_ca.csv'];%year,month,NC,SN,CC,SC
Tbl=readtable(nomefile);
[C,IA,IB] = intersect(years_study,years_williams);
BA=Tbl.NC+Tbl.SN;
for iyear=1:length(years_williams)
    i1=(iyear-1)*12+5
    i2=i1+4
    dum(iyear)=nansum(BA(i1:i2));
end
WILL_FOREST=zeros(length(years_study),1)*NaN;
WILL_FOREST(IA)=dum(IB); clear dum

%% Figure 1b


%MTBS 2021 Fire: Quarterly Release (August 10, 2022)
%MTBS has released 154 newly mapped fires from 2021 fire season for 16 CONUS states and Alaska.  It is the first quarterly release for 2021 season. With this data release to MTBS.gov, the total number of fires mapped by the project increases to 29,533.
%we skip 2021 data
[C,IA,IB] = intersect(years_study,years_frap);
FIRE=frap_forest_sierra_ncoast_year(IB)'/1000; 
col_fire=[0 0 0 ];
WILL_FOREST=WILL_FOREST/100000; col_will=[88 170 222]/255;
MTBS=MTBS*10; col_mtbs=[230 132 36]/255;
WEST=WEST/1000; col_west=[219 62 143]/255;

%
figure1 = figure; 
axes1 = axes('Parent',figure1);
plot(years_study,FIRE,'-o','color',col_fire,'LineWidth',1.5)
xlim([1970 2022])
hold on
plot(years_study,MTBS,'-s','color',col_mtbs,'LineWidth',1.5)
%plot(years_study,WILL_ALL,'-s','color',col_will,'LineWidth',1.5)
plot(years_study,WEST,'-d','color',col_west,'LineWidth',1.5)
plot(years_study,WILL_FOREST,'->','color',col_will,'LineWidth',1.5)
xlim([years_study(1)-1 years_study(end)+1])
legend('FRAP','MTBS','Westerling et al. (2003)','Williams et al. (2019)','Location','NorthWest')
plot(years_study,FIRE,'-o','color',col_fire,'LineWidth',1.5)
legend boxoff
%plot(years_study,FIRE,'-o','color',[255,0,0]/255,'LineWidth',1.5)
xlabel('Years','FontSize',20)
ylabel('Burned Area (1000 km^2)','FontSize',20)
set(gca,'FontSize',font_size)
file=[dir_out,'figure_S1.eps']
print( gcf, '-depsc2', file ,'-painters')

% 
[r p]=corrcoef(FIRE,MTBS,'rows','complete')
[r p]=corrcoef(FIRE,WEST,'rows','complete')
[r p]=corrcoef(FIRE,WILL_FOREST,'rows','complete')
