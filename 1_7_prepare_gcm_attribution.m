%%
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years_study=1971:2021;
% years_nat_all=1850:2020;
% years_nat=1950:2020;
years_scen=2015:2100;
years_hist=1950:2021;
base_period=1960:1982;
dir_gcm='~/Documents/dati/gcm_data_cal/'
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
best_start_tx=4;
best_stop_tx=10;
max_members=50;

%% load domain

filename = [dir_data, '/shapefiles/baileys/Baileys_ecoregions_cal_wgs84.shp'];
eco = shaperead(filename, 'UseGeoCoords', true); % imports lat/lon data
dum=eco.Lon';
xv=[dum;dum(1)]; clear dum
dum=eco.Lat';
yv=[dum;dum(1)]; clear dum

%% load obs

filename = [dir_data,'nclimgrid/nclimgrid.mat'];
load(filename);
TSMAX_obs = zeros(length(years_climate),1)*NaN;
for iyear=1:length(years_climate)
    i1 = (iyear - 1) * 12 + best_start_tx;
    i2 = (iyear - 1) * 12 + best_stop_tx;
    TSMAX_obs(iyear) = mean(nclimgrid(1,i1:i2));
end
%figure;plot(TSMAX_obs)
TSMAX_obs_anom=scale_base_period(TSMAX_obs,base_period,years_climate);
figure;plot(TSMAX_obs_anom)
%% hist-nat
%extract nome model and lon lat at 1??x1??
nomefile = ['/Users/marco/Documents/dati/gcm_data_cal/Turco-PNAS-2022-CMIP6/result_220119_01_tasmax-histall-histnat.mat'];
load(nomefile)
lat1=lat;
lon1=lon;
[LAT1, LON1] = meshgrid(lat1,lon1);     % new mesh
k=1;
for i=1:length(lon1)
    for j=1:length(lat1)
        coord(k,1)=lon1(i);
        coord(k,2)=lat1(j);
        k=k+1;
    end
end

% load hist-nat
nomefile = [dir_gcm,'result_220119_01_tasmax-histnat-12_all-member_nointerp_0121.mat'];
%nomefile = [dir_gcm,'result_220119_01_tasmax-histall-12_all-member_nointerp_0121.mat'];
load(nomefile)
%[member{12}]

% extract bias adjusted TSMAX
TSMAX_NAT = zeros([length(years_climate) length(model) max_members] )*NaN;            % preallocation
TSMAX_NAT_ANOM = zeros([length(years_climate) length(model) max_members] )*NaN;            % preallocation
kmod=0;
for imodel=1:length(model)
    model(imodel)
    time2=cell2mat(time(imodel));
    years_gcm=round(time2(1):time2(end)); 
    lat=cell2mat(latm(imodel));
    lon=cell2mat(lonm(imodel))-360;
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    [LAT, LON] = meshgrid(lat,lon);         % original mesh
    data_orig=cell2mat(tmax(imodel));
    [m,n,k,~] = size(data_orig);
    members=member{imodel};
    for imember=1:length(members)  
        members(imember)
        data_new = zeros([size(LAT1) k] )*NaN;            % preallocation
        for i = 1:k
            aux = data_orig(:,:,i,imember);
            data_new(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
        end 
        %Passo a matrix(month,space)
        tmp = NaN*zeros(size(data_new,3),size(data_new,1)*size(data_new,2));
        for im=1:size(data_new,3)
            ij=0;
            for i=1:size(data_new,1)
                for j=1:size(data_new,2)
                    ij=ij+1;
                    tmp(im,ij)=data_new(i,j,im);
                end
            end
        end
        % pr=nanmean(tmp,1)-273.15;
        % drawStations(coord,'size',1,'resolution','high','israster','true','color',pr','colormap',(jet))
        % %caxis([-3 3])
        % h=colorbar      
        % spatial average
        tx_nat(:)=nanmean(tmp(:,in),2)-273.15;%kelvin --> celsius
        %TSMAX=NaN*zeros(length(years_climate),1)
        for iyear=1:length(years_gcm)
            i1 = (iyear - 1) * 12 + best_start_tx;
            i2 = (iyear - 1) * 12 + best_stop_tx;
            dum(iyear) = mean(tx_nat(i1:i2));
        end
        clear tx_nat
        [C,IA,IB] = intersect(years_climate,years_gcm);
        TSMAX=NaN*zeros(length(years_climate),1);
        TSMAX(IA)=dum(IB);
        TSMAX(end)=nanmean(TSMAX(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        % bias adjustment
        [aux1,aux2]=bias_correction_ecdf...
            (TSMAX_obs,TSMAX,TSMAX);
        TSMAX_anom=scale_base_period(aux1,base_period,years_climate);
        %TSMAX_anom=aux1';
        TSMAX_NAT(:,imodel,imember)=aux1;
        TSMAX_NAT_ANOM(:,imodel,imember)=TSMAX_anom;
        kmod=kmod+1;
        model_nat(kmod)=strcat(model(imodel),';',members(imember));
        %figure;plot(TSMAX_anom); hold on; plot(TSMAX_obs_anom,'r')        
    end
end



TSMAX_NAT_ANOM_MIX = reshape(TSMAX_NAT_ANOM, [size(TSMAX_NAT_ANOM,1), size(TSMAX_NAT_ANOM,2)*size(TSMAX_NAT_ANOM,3)]);
TSMAX_NAT_MIX = reshape(TSMAX_NAT, [size(TSMAX_NAT,1), size(TSMAX_NAT,2)*size(TSMAX_NAT,3)]);
figure;plot(TSMAX_NAT_MIX,'r')
nnz(~isnan(TSMAX_NAT_MIX(1,:))) %--> 164 simulations
iok_nat=find(~isnan(TSMAX_NAT_MIX(1,:)));
clear TSMAX_NAT
TSMAX_NAT=TSMAX_NAT_MIX(:,iok_nat);
TSMAX_NAT_ANOM=TSMAX_NAT_ANOM_MIX(:,iok_nat);

%% hist-ssp245 1971-2021

%BCC-CSM2-MR_ssp245_r1i1p1f1
TSMAX_ALL = zeros([length(years_climate) size(TSMAX_NAT,2)] )*NaN;            % preallocation
TSMAX_ALL_ANOM = zeros([length(years_climate) size(TSMAX_NAT_ANOM,2)] )*NaN;            % preallocation
kmod=0;
kk=0;
for imodel=1:length(model)
    model(imodel)
    members=member{imodel};
    for imember=1:length(members)  
        kk=kk+1;
        members(imember)
        tx_hist = zeros([(length(years_hist)-7)*12] ,1)*NaN;            % preallocation
        tx_ssp245 = zeros([length(years_scen)*12],1 )*NaN;            % preallocation
        %tasmax_CMIP6Amon_INM-CM5-0_ssp585_r1i1p1f1_mon_2015-2100_CA_Yizhou.mat
        %namefile_hist=strcat(dir_gcm,'gcm_1_x_1/tasmax/','tasmax_CMIP6Amon_',model(imodel),'_historical_',members(imember),'_mon_1950-2014_CA_Yizhou.mat')
        %namefile_ssp245=strcat(dir_gcm,'gcm_1_x_1/tasmax/','tasmax_CMIP6Amon_',model(imodel),'_ssp245_',members(imember),'_mon_2015-2100_CA_Yizhou.mat')
        namefile_hist=strcat(dir_gcm,'gcm_1_x_1/tasmax/','tasmax_CMIP6Amon_',model(imodel),'_historical_',members(imember),'_mon.mat')
        namefile_ssp245=strcat(dir_gcm,'gcm_1_x_1/tasmax/','tasmax_CMIP6Amon_',model(imodel),'_ssp245_',members(imember),'_mon.mat')

        if (exist(char(namefile_hist), 'file') == 2 && exist(char(namefile_ssp245), 'file') == 2)
            load(char(namefile_hist))
            tx_hist=tx_gcm;
            clear tx_gcm
            load(char(namefile_ssp245))
            tx_ssp245=tx_gcm;
            if (length(tx_ssp245)>=84)
                kmod=kmod+1;
                model_com(kmod)=strcat(model(imodel),';',members(imember))
                model_member_num(kmod,1)=imodel;
                model_member_num(kmod,2)=imember;
                clear tx_gcm
            else
                continue
            end
        else
            continue
        end        
        tx_gcm=[tx_hist; tx_ssp245(1:84)]; %hist up to 2014, ssp245 up to 2021, 84 months
        for iyear=1:length(years_hist)
            i1 = (iyear - 1) * 12 + best_start_tx;
            i2 = (iyear - 1) * 12 + best_stop_tx;
            TSMAX(iyear) = mean(tx_gcm(i1:i2));
        end
        %figure;plot(tx_gcm)
        %figure;plot(TSMAX)
        clear tx_gcm
        [C,IA,IB] = intersect(years_climate,years_hist);
        TSMAX=TSMAX(IB);
        % bias adjustment
        [aux1,aux2]=bias_correction_ecdf...
            (TSMAX_obs,TSMAX,TSMAX);
        TSMAX_anom=scale_base_period(aux1,base_period,years_climate);
        %TSMAX_anom=aux1';
        TSMAX_ALL(IA,kk)=aux1';
        TSMAX_ALL_ANOM(IA,kk)=TSMAX_anom;
    end
end

figure;plot(TSMAX_ALL,'r')
nnz(~isnan(TSMAX_ALL(1,:))) %--> 86 simulations
iok_all=find(~isnan(TSMAX_ALL(1,:)));

%% save outputs
TSMAX_ALL=TSMAX_ALL(:,iok_all);
TSMAX_NAT=TSMAX_NAT(:,iok_all);
TSMAX_ALL_ANOM=TSMAX_ALL_ANOM(:,iok_all);
TSMAX_NAT_ANOM=TSMAX_NAT_ANOM(:,iok_all);
model_com
length(model_com(:))
model_member_num
filename = [dir_data 'gcms/GCMs_attribution.mat'];
save(filename,'model_com','model_member_num','TSMAX_ALL','TSMAX_NAT','TSMAX_ALL_ANOM','TSMAX_NAT_ANOM')  

%% prec

%%
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years_study=1971:2021;
% years_nat_all=1850:2020;
% years_nat=1950:2020;
years_scen=2015:2100;
years_hist=1950:2021;
base_period=1960:1982;
dir_gcm='~/Documents/dati/gcm_data_cal/'
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
best_start_pr=4;
best_stop_pr=10;
max_members=50;

%% load domain

filename = [dir_data, '/shapefiles/baileys/Baileys_ecoregions_cal_wgs84.shp'];
eco = shaperead(filename, 'UseGeoCoords', true); % imports lat/lon data
dum=eco.Lon';
xv=[dum;dum(1)]; clear dum
dum=eco.Lat';
yv=[dum;dum(1)]; clear dum

%% load obs

filename = [dir_data,'nclimgrid/nclimgrid.mat'];
load(filename);

PREC_obs = zeros(length(years_climate),1)*NaN;
for iyear=1:length(years_climate)
    i1 = (iyear - 1) * 12 + best_start_pr;
    i2 = (iyear - 1) * 12 + best_stop_pr;
    PREC_obs(iyear) = sum(nclimgrid(2,i1:i2));
end
%figure;plot(TSMAX_obs)
%PREC_obs=scale_base_period(PREC_obs,base_period,years_climate);
figure;plot(PREC_obs)

%% hist-nat
%extract nome model and lon lat at 1??x1??
nomefile = ['/Users/marco/Documents/dati/gcm_data_cal/Turco-PNAS-2022-CMIP6/result_220119_01_tasmax-histall-histnat.mat'];
load(nomefile)
lat1=lat;
lon1=lon;
[LAT1, LON1] = meshgrid(lat1,lon1);     % new mesh
k=1;
for i=1:length(lon1)
    for j=1:length(lat1)
        coord(k,1)=lon1(i);
        coord(k,2)=lat1(j);
        k=k+1;
    end
end

% load hist-nat
nomefile = ['/Users/marco/Documents/dati/gcm_data_cal/Turco-PNAS-2022-CMIP6/result_221203_01_pr-histnat_12_nointerp.mat'];
load(nomefile)
%[member{12}]

% extract bias adjusted TSMAX
PREC_NAT = zeros([length(years_climate) length(models_u) max_members] )*NaN;            % preallocation
PREC_NAT_ANOM = zeros([length(years_climate) length(model) max_members] )*NaN;            % preallocation
kmod=0;
for imodel=1:length(model)
    models_u(imodel)
    time2=cell2mat(time{1,imodel});
    years_gcm=round(time2(1):time2(end)); 
    lat=cell2mat(latm(imodel));
    lon=cell2mat(lonm(imodel))-360;
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    [LAT, LON] = meshgrid(lat,lon);         % original mesh
    members=variants_u{imodel};
    for imember=1:length(members)  
        members(imember)
        data_orig=(Pr{1,imodel}{1,imember});
        [m,n,k] = size(data_orig);
        data_new = zeros([size(LAT1) k] )*NaN;            % preallocation
        for i = 1:k
            aux = data_orig(:,:,i);
            data_new(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
        end 
        %Passo a matrix(month,space)
        tmp = NaN*zeros(size(data_new,3),size(data_new,1)*size(data_new,2));
        for im=1:size(data_new,3)
            ij=0;
            for i=1:size(data_new,1)
                for j=1:size(data_new,2)
                    ij=ij+1;
                    tmp(im,ij)=data_new(i,j,im);
                end
            end
        end
        % pr=nanmean(tmp,1)-273.15;
        % drawStations(coord,'size',1,'resolution','high','israster','true','color',pr','colormap',(jet))
        % %caxis([-3 3])
        % h=colorbar      
        % spatial average
        pr_nat(:)=nanmean(tmp(:,in),2);
        %TSMAX=NaN*zeros(length(years_climate),1)
        for iyear=1:length(years_gcm)
            i1 = (iyear - 1) * 12 + best_start_pr;
            i2 = (iyear - 1) * 12 + best_stop_pr;
            dum(iyear) = mean(pr_nat(i1:i2));
        end
        clear pr_nat
        [C,IA,IB] = intersect(years_climate,years_gcm);
        PREC=NaN*zeros(length(years_climate),1);
        PREC(IA)=dum(IB);
        PREC(end)=nanmean(PREC(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        % bias adjustment
        [aux1,aux2]=bias_correction_ecdf...
            (PREC_obs,PREC,PREC);
        PREC_anom=scale_base_period(aux1,base_period,years_climate);
        %TSMAX_anom=aux1';
        PREC_NAT(:,imodel,imember)=aux1;
        PREC_NAT_ANOM(:,imodel,imember)=PREC_anom;
        kmod=kmod+1;
        model_nat(kmod)=strcat(model(imodel),';',members(imember));
        %figure;plot(TSMAX_anom); hold on; plot(TSMAX_obs_anom,'r')        
    end
end

PREC_NAT_ANOM_MIX = reshape(PREC_NAT_ANOM, [size(PREC_NAT_ANOM,1), size(PREC_NAT_ANOM,2)*size(PREC_NAT_ANOM,3)]);
PREC_NAT_MIX = reshape(PREC_NAT, [size(PREC_NAT,1), size(PREC_NAT,2)*size(PREC_NAT,3)]);
figure;plot(PREC_NAT_MIX,'r')
nnz(~isnan(PREC_NAT_MIX(1,:))) %--> 164 simulations
iok_nat=find(~isnan(PREC_NAT_MIX(1,:)));
clear TSMAX_NAT
PREC_NAT=PREC_NAT_MIX(:,iok_nat);
PREC_NAT_ANOM=PREC_NAT_ANOM_MIX(:,iok_nat);

%% hist-ssp245 1971-2021

%BCC-CSM2-MR_ssp245_r1i1p1f1
PREC_ALL = zeros([length(years_climate) size(PREC_NAT,2)] )*NaN;            % preallocation
PREC_ALL_ANOM = zeros([length(years_climate) size(PREC_NAT_ANOM,2)] )*NaN;            % preallocation
kmod=0;
kk=0;
for imodel=1:length(models_u)
    models_u(imodel)
    members=variants_u{imodel};
    for imember=1:length(members)  
        kk=kk+1;
        members(imember)
        pr_hist = zeros([(length(years_hist)-7)*12] ,1)*NaN;            % preallocation
        pr_ssp245 = zeros([length(years_scen)*12],1 )*NaN;            % preallocation
        %tasmax_CMIP6Amon_INM-CM5-0_ssp585_r1i1p1f1_mon_2015-2100_CA_Yizhou.mat
        %namefile_hist=strcat(dir_gcm,'gcm_1_x_1/pr/','pr_CMIP6Amon_',models_u(imodel),'_historical_',members(imember),'_mon_1950-2014_CA_Yizhou.mat')
        %namefile_ssp245=strcat(dir_gcm,'gcm_1_x_1/pr/','pr_CMIP6Amon_',models_u(imodel),'_ssp245_',members(imember),'_mon_2015-2100_CA_Yizhou.mat')
        namefile_hist=strcat(dir_gcm,'gcm_1_x_1/pr/','pr_CMIP6Amon_',models_u(imodel),'_historical_',members(imember),'_mon.mat')
        namefile_ssp245=strcat(dir_gcm,'gcm_1_x_1/pr/','pr_CMIP6Amon_',models_u(imodel),'_ssp245_',members(imember),'_mon.mat')
        if (exist(char(namefile_hist), 'file') == 2 && exist(char(namefile_ssp245), 'file') == 2)
            load(char(namefile_hist))
            pr_hist=pr_gcm;
            clear pr_gcm
            load(char(namefile_ssp245))
            pr_ssp245=pr_gcm;
            if (length(pr_ssp245)>=84)
                kmod=kmod+1;
                model_com(kmod)=strcat(model(imodel),';',members(imember))
                model_member_num(kmod,1)=imodel;
                model_member_num(kmod,2)=imember;
                clear pr_gcm
            else
                continue
            end
        else
            continue
        end        
        pr_gcm=[pr_hist; pr_ssp245(1:84)]; %hist up to 2014, ssp245 up to 2021, 84 months
        for iyear=1:length(years_hist)
            i1 = (iyear - 1) * 12 + best_start_pr;
            i2 = (iyear - 1) * 12 + best_stop_pr;
            PREC(iyear) = mean(pr_gcm(i1:i2));
        end
        %figure;plot(tx_gcm)
        %figure;plot(TSMAX)
        clear pr_gcm
        [C,IA,IB] = intersect(years_climate,years_hist);
        PREC=PREC(IB);
        % bias adjustment
        [aux1,aux2]=bias_correction_ecdf...
            (PREC_obs,PREC,PREC);
        PREC_anom=scale_base_period(aux1,base_period,years_climate);
        %TSMAX_anom=aux1';
        PREC_ALL(IA,kk)=aux1';
        PREC_ALL_ANOM(IA,kk)=PREC_anom;
    end
end

figure;plot(PREC_obs)
figure;plot(PREC_ALL,'r')
nnz(~isnan(PREC_ALL(1,:))) %--> 86 simulations
iok_all=find(~isnan(PREC_ALL(1,:)));


%% save outputs
PREC_ALL=PREC_ALL(:,iok_all);
PREC_NAT=PREC_NAT(:,iok_all);
PREC_ALL_ANOM=PREC_ALL_ANOM(:,iok_all);
PREC_NAT_ANOM=PREC_NAT_ANOM(:,iok_all);

figure;plot(PREC_ALL_ANOM,'r')
figure;plot(PREC_ALL,'r')

model_com(:)
model_member_num
filename = [dir_data 'gcms/GCMs_attribution_prec.mat'];
save(filename,'model_com','model_member_num','PREC_ALL','PREC_NAT','PREC_ALL_ANOM','PREC_NAT_ANOM')  

%% VPD
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years_study=1971:2021;
% years_nat_all=1850:2020;
% years_nat=1950:2020;
years_scen=2015:2100;
years_hist=1950:2021;
base_period=1960:1982;

dir_gcm='~/Documents/dati/gcm_data_cal/Turco-PNAS-2022-CMIP6/'
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
best_start_tx=4;
best_stop_tx=10;
max_members=50;

%% load domain
filename = [dir_data, '/shapefiles/baileys/Baileys_ecoregions_cal_wgs84.shp'];
eco = shaperead(filename, 'UseGeoCoords', true); % imports lat/lon data
dum=eco.Lon';
xv=[dum;dum(1)]; clear dum
dum=eco.Lat';
yv=[dum;dum(1)]; clear dum

%% load obs
%climatic variables prism
filename = [dir_data,'prism/vpd_prism.mat']; %period 1960:2021
load(filename);

VPD_obs = zeros(length(years_climate),1)*NaN;
for iyear=1:length(years_climate)
    i1 = (iyear - 1) * 12 + best_start_tx;
    i2 = (iyear - 1) * 12 + best_stop_tx;
    VPD_obs(iyear) = mean(vpd_prism(i1:i2));
end
VPD_obs_anom=scale_base_period(VPD_obs,base_period,years_climate);
%% figure;plot(VPD_obs)

%% hist-nat
%extract nome model and lon lat at 1??x1??
nomefile = [dir_gcm,'result_220119_01_tasmax-histall-histnat.mat'];
load(nomefile)
lat1=lat;
lon1=lon;
[LAT1, LON1] = meshgrid(lat1,lon1);     % new mesh
auxLat=LAT1';auxLon=LON1';
coord=[auxLon(:) auxLat(:)];

%% nameFiles={'result_221203_01_rhmean-histnat_12_nointerp.mat';'result_221203_01_tmax-histnat_12_nointerp.mat';'result_221203_01_tmin-histnat_12_nointerp.mat';...
%% 'result_221203_01_rhmean-hist_34_nointerp.mat';'result_221203_01_tmax-hist_34_nointerp.mat';'result_221203_01_tmin-hist_34_nointerp.mat';...
%% 'result_221212_01_rhmean-ssp245_28_nointerp.mat';'result_221212_01_tmax-ssp245_28_nointerp.mat';'result_221212_01_tmin-ssp245_28_nointerp.mat'};
%% for i=1:3
%% 	nomefile = [dir_gcm nameFiles{i}];load(nomefile)
%% 	nomefile = [dir_gcm nameFiles{3+i}];aux=load(nomefile);
%% 	models_u=intersect(deblank(models_u),deblank(aux.models_u));
%% 	nomefile = [dir_gcm nameFiles{2*3+i}];aux=load(nomefile);
%% 	models_u=intersect(deblank(models_u),deblank(aux.models_u));
%% 	disp(nameFiles{i})
%% 	strvcat(models_u)
%% end
gcmNames={'ACCESS-CM2';'ACCESS-ESM1-5';'CNRM-CM6-1';'CanESM5';'FGOALS-g3';'GFDL-ESM4';'GISS-E2-1-G';'HadGEM3-GC31-LL';'IPSL-CM6A-LR';'MIROC6';'MRI-ESM2-0'};

%% load hist-nat
nomefile = [dir_gcm,'result_221203_01_rhmean-histnat_12_nointerp.mat'];
load(nomefile);
[models_u,I1,I2]=intersect(deblank(models_u),deblank(gcmNames));
Rhmean=Rhmean(I1);latm=latm(I1);lonm=lonm(I1);
time=time(I1);variants_u=variants_u(I1);
models_u_rhmean_histnat=models_u;
time_rhmean_histnat=time;
variants_u_rhmean_histnat=variants_u;

%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221203_01_tmax-histnat_12_nointerp.mat'];
nomefile = [dir_gcm,'/result_221203_01_tmax-histnat_12_nointerp.mat'];
load(nomefile)
Tmax=Tmax(I1);latm=latm(I1);lonm=lonm(I1);
time=time(I1);variants_u=variants_u(I1);models_u=models_u(I1);
models_u_tmax_histnat=models_u;
time_tmax_histnat=time;
variants_u_tmax_histnat=variants_u;

%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221203_01_tmin-histnat_12_nointerp.mat'];
nomefile = [dir_gcm,'/result_221203_01_tmin-histnat_12_nointerp.mat'];
load(nomefile)
Tmin=Tmin(I1);latm=latm(I1);lonm=lonm(I1);
time=time(I1);variants_u=variants_u(I1);models_u=models_u(I1);
models_u_tmin_histnat=models_u;
time_tmin_histnat=time;
variants_u_tmin_histnat=variants_u;

%celldisp(variants_u_rhmean_histnat)
%[member{12}]
% extract bias adjusted TSMAX
VPD_NAT = zeros([length(years_climate) length(models_u) max_members] )*NaN;            % preallocation
VPD_NAT_ANOM = zeros([length(years_climate) length(model) max_members] )*NaN;            % preallocation
kmod=0;
for imodel=1:length(models_u_rhmean_histnat)
    disp(models_u_rhmean_histnat{imodel})
    time2=cell2mat(time_rhmean_histnat{1,imodel});
    years_gcm=round(time2(1):time2(end)); 
    lat=cell2mat(latm(imodel));
    lon=cell2mat(lonm(imodel))-360;
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    [LAT, LON] = meshgrid(lat,lon);         % original mesh
    members=variants_u_rhmean_histnat{imodel};
    for imember=1:length(members)  
        %% members(imember);
        data_orig_rhmean=(Rhmean{1,imodel}{1,imember});
        data_orig_tmax=(Tmax{1,imodel}{1,imember});
        data_orig_tmin=(Tmin{1,imodel}{1,imember});
        [m,n,k] = size(data_orig_rhmean);
        data_new_rhmean = zeros([size(LAT1) k] )*NaN;            % preallocation
        data_new_tmax = zeros([size(LAT1) k] )*NaN;            % preallocation
        data_new_tmin = zeros([size(LAT1) k] )*NaN;            % preallocation
        for i = 1:k
            aux = data_orig_rhmean(:,:,i);
            data_new_rhmean(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
            aux = data_orig_tmax(:,:,i);
            data_new_tmax(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
            aux = data_orig_tmin(:,:,i);
            data_new_tmin(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
        end 
        %Passo a matrix(month,space)
        tmp_rhmean = NaN*zeros(size(data_new_rhmean,3),size(data_new_rhmean,1)*size(data_new_rhmean,2));
        tmp_tmax = NaN*zeros(size(data_new_rhmean,3),size(data_new_rhmean,1)*size(data_new_rhmean,2));
        tmp_tmin = NaN*zeros(size(data_new_rhmean,3),size(data_new_rhmean,1)*size(data_new_rhmean,2));
        for im=1:size(data_new_rhmean,3)
            ij=0;
            for i=1:size(data_new_rhmean,1)
                for j=1:size(data_new_rhmean,2)
                    ij=ij+1;
                    tmp_rhmean(im,ij)=data_new_rhmean(i,j,im);
                    tmp_tmax(im,ij)=data_new_tmax(i,j,im);
                    tmp_tmin(im,ij)=data_new_tmin(i,j,im);
                end
            end
        end
        % spatial average
        rhmean_nat(:)=nanmean(tmp_rhmean(:,in),2);
        tmax_nat(:)=nanmean(tmp_tmax(:,in),2);
        tmin_nat(:)=nanmean(tmp_tmin(:,in),2);
        for iyear=1:length(years_gcm)
            i1 = (iyear - 1) * 12 + best_start_tx;
            i2 = (iyear - 1) * 12 + best_stop_tx;
            dum_rhmean(iyear) = mean(rhmean_nat(i1:i2));
            dum_tmax(iyear) = mean(tmax_nat(i1:i2));
            dum_tmin(iyear) = mean(tmin_nat(i1:i2));
        end
        clear rhmean_nat tmax_nat tmin_nat
        [C,IA,IB] = intersect(years_climate,years_gcm);
        RHMEAN=NaN*zeros(length(years_climate),1);
        RHMEAN(IA)=dum_rhmean(IB);
        RHMEAN(end)=nanmean(RHMEAN(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        
        TMAX=NaN*zeros(length(years_climate),1);
        TMAX(IA)=dum_tmax(IB);
        TMAX(end)=nanmean(TMAX(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        TMAX=TMAX-273.15;
        
        TMIN=NaN*zeros(length(years_climate),1);
        TMIN(IA)=dum_tmin(IB);
        TMIN(end)=nanmean(TMIN(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        TMIN=TMIN-273.15;
        
        %esmn=0.6108 * exp((17.27 * TMIN) ./ (TMIN + 237.3))
        %calculate VPD
        VPD=calculate_vpd_gcms(RHMEAN,TMAX,TMIN);
        %plot(VPD)
        
        % bias adjustment
        [aux1,aux2]=bias_correction_ecdf(VPD_obs,VPD,VPD);
        VPD_anom=scale_base_period(aux1,base_period,years_climate);
        %TSMAX_anom=aux1';
        VPD_NAT(:,imodel,imember)=aux1;
        VPD_NAT_ANOM(:,imodel,imember)=VPD_anom;
        kmod=kmod+1;
        model_nat(kmod)=strcat(models_u_rhmean_histnat(imodel),';',members(imember));
        %figure;plot(TSMAX_anom); hold on; plot(TSMAX_obs_anom,'r')        
    end
end
VPD_NAT_ANOM_MIX = reshape(VPD_NAT_ANOM, [size(VPD_NAT_ANOM,1), size(VPD_NAT_ANOM,2)*size(VPD_NAT_ANOM,3)]);
VPD_NAT_MIX = reshape(VPD_NAT, [size(VPD_NAT,1), size(VPD_NAT,2)*size(VPD_NAT,3)]);
%% figure;plot(VPD_NAT_MIX,'r')
nnz(~isnan(VPD_NAT_MIX(1,:))) %--> 164 simulations
iok_nat=find(~isnan(VPD_NAT_MIX(1,:)));
clear VPD_NAT
VPD_NAT=VPD_NAT_MIX(:,iok_nat);
VPD_NAT_ANOM=VPD_NAT_ANOM_MIX(:,iok_nat);

%% hist-ssp245 1971-2021
%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221203_01_rhmean-hist_34_nointerp.mat'];
nomefile = [dir_gcm,'/result_221203_01_rhmean-hist_34_nointerp.mat'];
auxh=load(nomefile);
[models_u,I1,I2]=intersect(deblank(auxh.models_u),deblank(gcmNames));
auxh.Rhmean=auxh.Rhmean(I1);
auxh.latm=auxh.latm(I1);
auxh.lonm=auxh.lonm(I1);
auxh.models_u=auxh.models_u(I1);
auxh.time=auxh.time(I1);
auxh.variants_u=auxh.variants_u(I1);

nomefile = [dir_gcm,'/result_221212_01_rhmean-ssp245_28_nointerp.mat'];
aux=load(nomefile);
[models_u,I1,I2]=intersect(deblank(aux.models_u),deblank(gcmNames));
aux.Rhmean=aux.Rhmean(I1);
aux.latm=aux.latm(I1);
aux.lonm=aux.lonm(I1);
aux.models_u=aux.models_u(I1);
aux.time=aux.time(I1);
aux.variants_u=aux.variants_u(I1);

Rhmean=cell(1,length(aux.Rhmean));
variants_u=cell(1,length(aux.Rhmean));
time=cell(1,length(aux.Rhmean));
latm=cell(1,length(aux.Rhmean));
lonm=cell(1,length(aux.Rhmean));
for i=1:length(aux.Rhmean)
	[auxMembers,J1,J2]=intersect(deblank(auxh.variants_u{i}),deblank(aux.variants_u{i}));
	Rhmean{i}=cell(1,length(auxMembers));
	time{i}=cell(1,length(auxMembers));
	variants_u{i}=auxMembers;
	latm{i}=auxh.latm{i};
	lonm{i}=auxh.lonm{i};
	for j=1:length(auxMembers)
		time{i}{j}=[auxh.time{i}{J1(j)};aux.time{i}{J2(j)}];
		data=repmat(NaN,[size(aux.Rhmean{i}{J2(j)},1) size(aux.Rhmean{i}{J2(j)},2) length(time{i}{j})]);
		data(:,:,1:length(auxh.time{i}{J1(j)}))=auxh.Rhmean{i}{J2(j)};
		data(:,:,(1+length(auxh.time{i}{J1(j)})):end)=aux.Rhmean{i}{J2(j)};
		Rhmean{i}{j}=data;clear data
	end
end
models_u_rhmean=models_u;
time_rhmean=time;
variants_u_rhmean=variants_u;

%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221203_01_tmin-hist_34_nointerp.mat'];
%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221212_01_tmin-ssp245_28_nointerp.mat'];
nomefile = [dir_gcm,'/result_221203_01_tmin-hist_34_nointerp.mat'];
auxh=load(nomefile);
[models_u,I1,I2]=intersect(deblank(auxh.models_u),deblank(gcmNames));
auxh.Tmin=auxh.Tmin(I1);
auxh.latm=auxh.latm(I1);
auxh.lonm=auxh.lonm(I1);
auxh.models_u=auxh.models_u(I1);
auxh.time=auxh.time(I1);
auxh.variants_u=auxh.variants_u(I1);

nomefile = [dir_gcm,'/result_221212_01_tmin-ssp245_28_nointerp.mat'];
aux=load(nomefile);
[models_u,I1,I2]=intersect(deblank(aux.models_u),deblank(gcmNames));
aux.Tmin=aux.Tmin(I1);
aux.latm=aux.latm(I1);
aux.lonm=aux.lonm(I1);
aux.models_u=aux.models_u(I1);
aux.time=aux.time(I1);
aux.variants_u=aux.variants_u(I1);

Tmin=cell(1,length(aux.Tmin));
variants_u=cell(1,length(aux.Tmin));
time=cell(1,length(aux.Tmin));
latm=cell(1,length(aux.Tmin));
lonm=cell(1,length(aux.Tmin));
for i=1:length(aux.Tmin)
	[auxMembers,J1,J2]=intersect(deblank(auxh.variants_u{i}),deblank(aux.variants_u{i}));
	Tmin{i}=cell(1,length(auxMembers));
	time{i}=cell(1,length(auxMembers));
	variants_u{i}=auxMembers;
	latm{i}=auxh.latm{i};
	lonm{i}=auxh.lonm{i};
	for j=1:length(auxMembers)
		time{i}{j}=[auxh.time{i}{J1(j)};aux.time{i}{J2(j)}];
		data=repmat(NaN,[size(aux.Tmin{i}{J2(j)},1) size(aux.Tmin{i}{J2(j)},2) length(time{i}{j})]);
		data(:,:,1:length(auxh.time{i}{J1(j)}))=auxh.Tmin{i}{J2(j)};
		data(:,:,(1+length(auxh.time{i}{J1(j)})):end)=aux.Tmin{i}{J2(j)};
		Tmin{i}{j}=data;clear data
	end
end
models_u_tmin=models_u;
time_tmin=time;
variants_u_tmin=variants_u;

%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221203_01_tmax-hist_34_nointerp.mat'];
%% nomefile = [dir_gcm,'Turco-PNAS-2022-CMIP6/result_221212_01_tmax-ssp245_28_nointerp.mat'];
nomefile = [dir_gcm,'/result_221203_01_tmax-hist_34_nointerp.mat'];
auxh=load(nomefile);
[models_u,I1,I2]=intersect(deblank(auxh.models_u),deblank(gcmNames));
auxh.Tmax=auxh.Tmax(I1);
auxh.latm=auxh.latm(I1);
auxh.lonm=auxh.lonm(I1);
auxh.models_u=auxh.models_u(I1);
auxh.time=auxh.time(I1);
auxh.variants_u=auxh.variants_u(I1);

nomefile = [dir_gcm,'/result_221212_01_tmax-ssp245_28_nointerp.mat'];
aux=load(nomefile);
[models_u,I1,I2]=intersect(deblank(aux.models_u),deblank(gcmNames));
aux.Tmax=aux.Tmax(I1);
aux.latm=aux.latm(I1);
aux.lonm=aux.lonm(I1);
aux.models_u=aux.models_u(I1);
aux.time=aux.time(I1);
aux.variants_u=aux.variants_u(I1);

Tmax=cell(1,length(aux.Tmax));
variants_u=cell(1,length(aux.Tmax));
time=cell(1,length(aux.Tmax));
latm=cell(1,length(aux.Tmax));
lonm=cell(1,length(aux.Tmax));
for i=1:length(aux.Tmax)
	[auxMembers,J1,J2]=intersect(deblank(auxh.variants_u{i}),deblank(aux.variants_u{i}));
	Tmax{i}=cell(1,length(auxMembers));
	time{i}=cell(1,length(auxMembers));
	variants_u{i}=auxMembers;
	latm{i}=auxh.latm{i};
	lonm{i}=auxh.lonm{i};
	for j=1:length(auxMembers)
		time{i}{j}=[auxh.time{i}{J1(j)};aux.time{i}{J2(j)}];
		data=repmat(NaN,[size(aux.Tmax{i}{J2(j)},1) size(aux.Tmax{i}{J2(j)},2) length(time{i}{j})]);
		data(:,:,1:length(auxh.time{i}{J1(j)}))=auxh.Tmax{i}{J2(j)};
		data(:,:,(1+length(auxh.time{i}{J1(j)})):end)=aux.Tmax{i}{J2(j)};
		Tmax{i}{j}=data;clear data
	end
end
models_u_tmax=models_u;
time_tmax=time;
variants_u_tmax=variants_u;

%[C,ia,ib_hist] = intersect(models_u_rhmean_histnat,models_u);
%[C,ia,ib_ssp245] = intersect(models_u_rhmean_histnat,models_u);
%[models_u_rhmean_ssp245' models_u_tmin_ssp245' models_u_tmax_ssp245']
%[models_u_rhmean_hist' models_u_tmin_hist' models_u_tmax_hist']

% extract bias adjusted TSMAX
VPD_ALL = zeros([length(years_climate) length(models_u) max_members] )*NaN;            % preallocation
VPD_ALL_ANOM = zeros([length(years_climate) length(models_u) max_members] )*NaN;            % preallocation
kmod=0;
for imodel=1:length(models_u_rhmean)
    disp(models_u_rhmean{imodel})
	if imodel==7,
		time2=cell2mat(time_rhmean{1,imodel}(2:end));
		members=intersect(intersect(variants_u_rhmean{imodel}(2:end),variants_u_tmin{imodel}(2:end)),variants_u_tmax{imodel}(2:end));
	elseif imodel==10
		indTime=zeros(length(variants_u_rhmean{imodel}),1);for iTime=1:length(variants_u_rhmean{imodel}),if length(time_rhmean{1,imodel}{iTime})>=3000,indTime(iTime)=1;end;end
		time2=cell2mat(time_rhmean{1,imodel}(find(indTime==1)));
		members=intersect(intersect(variants_u_rhmean{imodel}(find(indTime==1)),variants_u_tmin{imodel}),variants_u_tmax{imodel});
	else
		time2=cell2mat(time_rhmean{1,imodel});
		members=intersect(intersect(variants_u_rhmean{imodel},variants_u_tmin{imodel}),variants_u_tmax{imodel});
	end
    years_gcm=round(time2(1):time2(end)); 
    lat=cell2mat(latm(imodel));
    lon=cell2mat(lonm(imodel))-360;
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    [LAT, LON] = meshgrid(lat,lon);         % original mesh
    [auxM,M1,MR]=intersect(deblank(members),deblank(variants_u_rhmean{imodel}));
    [auxM,M1,MN]=intersect(deblank(members),deblank(variants_u_tmin{imodel}));
    [auxM,M1,MX]=intersect(deblank(members),deblank(variants_u_tmax{imodel}));
    for imember=1:length(members)  
        %% members(imember);
        data_orig_rhmean=(Rhmean{1,imodel}{1,MR(imember)});
        data_orig_tmax=(Tmax{1,imodel}{1,MX(imember)});
        data_orig_tmin=(Tmin{1,imodel}{1,MN(imember)});
        [m,n,k] = size(data_orig_rhmean);
        data_new_rhmean = zeros([size(LAT1) k] )*NaN;            % preallocation
        data_new_tmax = zeros([size(LAT1) k] )*NaN;            % preallocation
        data_new_tmin = zeros([size(LAT1) k] )*NaN;            % preallocation
        for i = 1:k
            aux = data_orig_rhmean(:,:,i);
            data_new_rhmean(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
            aux = data_orig_tmax(:,:,i);
            data_new_tmax(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
            aux = data_orig_tmin(:,:,i);
            data_new_tmin(:,:,i) = interp2(LAT,LON,aux,LAT1,LON1);
        end 
        %Passo a matrix(month,space)
        tmp_rhmean = NaN*zeros(size(data_new_rhmean,3),size(data_new_rhmean,1)*size(data_new_rhmean,2));
        tmp_tmax = NaN*zeros(size(data_new_rhmean,3),size(data_new_rhmean,1)*size(data_new_rhmean,2));
        tmp_tmin = NaN*zeros(size(data_new_rhmean,3),size(data_new_rhmean,1)*size(data_new_rhmean,2));
        for im=1:size(data_new_rhmean,3)
            ij=0;
            for i=1:size(data_new_rhmean,1)
                for j=1:size(data_new_rhmean,2)
                    ij=ij+1;
                    tmp_rhmean(im,ij)=data_new_rhmean(i,j,im);
                    tmp_tmax(im,ij)=data_new_tmax(i,j,im);
                    tmp_tmin(im,ij)=data_new_tmin(i,j,im);
                end
            end
        end
        % spatial average
        rhmean_nat(:)=nanmean(tmp_rhmean(:,in),2);
        tmax_nat(:)=nanmean(tmp_tmax(:,in),2);
        tmin_nat(:)=nanmean(tmp_tmin(:,in),2);
        for iyear=1:length(years_gcm)
            i1 = (iyear - 1) * 12 + best_start_tx;
            i2 = (iyear - 1) * 12 + best_stop_tx;
            dum_rhmean(iyear) = mean(rhmean_nat(i1:i2));
            dum_tmax(iyear) = mean(tmax_nat(i1:i2));
            dum_tmin(iyear) = mean(tmin_nat(i1:i2));
        end
        clear rhmean_nat tmax_nat tmin_nat
        [C,IA,IB] = intersect(years_climate,years_gcm);
        RHMEAN=NaN*zeros(length(years_climate),1);
        RHMEAN(IA)=dum_rhmean(IB);
        RHMEAN(end)=nanmean(RHMEAN(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        
        TMAX=NaN*zeros(length(years_climate),1);
        TMAX(IA)=dum_tmax(IB);
        TMAX(end)=nanmean(TMAX(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        TMAX=TMAX-273.15;
        
        TMIN=NaN*zeros(length(years_climate),1);
        TMIN(IA)=dum_tmin(IB);
        TMIN(end)=nanmean(TMIN(42:61)); %for 2021 i used the mean of 2001-2020 data of hist-nat runs
        TMIN=TMIN-273.15;
        
        %esmn=0.6108 * exp((17.27 * TMIN) ./ (TMIN + 237.3))
        %calculate VPD
        VPD=calculate_vpd_gcms(RHMEAN,TMAX,TMIN);
        %plot(VPD)
        
        % bias adjustment
        [aux1,aux2]=bias_correction_ecdf(VPD_obs,VPD,VPD);
        VPD_anom=scale_base_period(aux1,base_period,years_climate);
        %TSMAX_anom=aux1';
        VPD_ALL(:,imodel,imember)=aux1;
        VPD_ALL_ANOM(:,imodel,imember)=VPD_anom;
        kmod=kmod+1;
        model_all(kmod)=strcat(models_u_rhmean(imodel),';',members(imember));
        %figure;plot(TSMAX_anom); hold on; plot(TSMAX_obs_anom,'r')        
    end
end
VPD_ALL_ANOM_MIX = reshape(VPD_ALL_ANOM, [size(VPD_ALL_ANOM,1), size(VPD_ALL_ANOM,2)*size(VPD_ALL_ANOM,3)]);
VPD_ALL_MIX = reshape(VPD_ALL, [size(VPD_ALL,1), size(VPD_ALL,2)*size(VPD_ALL,3)]);
%% figure;plot(VPD_NAT_MIX,'r')
nnz(~isnan(VPD_ALL_MIX(1,:))) %--> 164 simulations
iok_all=find(~isnan(VPD_ALL_MIX(1,:)));
clear VPD_ALL
VPD_ALL=VPD_ALL_MIX(:,iok_all);
VPD_ALL_ANOM=VPD_ALL_ANOM_MIX(:,iok_all);

%% save outputs
filename = [dir_data 'gcms/GCMs_attribution_vpd.mat'];
save(filename,'model_nat','model_all','VPD_ALL','VPD_NAT','VPD_ALL_ANOM','VPD_NAT_ANOM')  
