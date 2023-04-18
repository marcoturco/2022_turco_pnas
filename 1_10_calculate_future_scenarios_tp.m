%% plot figure 2b, 2c, 2d
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_study=1971:2021;
base_period=1995:2014;
years_frap=1950:2021;
years_hist = 1950:2014; %years hist-nat simulations
years_fut = 2015:2100; %years future simulations
years_sim=1950:2100;
years_sim_final=1971:2100;
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
alpha = 0.05;
NB=10000;
best_start_tx=4;
best_stop_tx=10;
font_size=12;

%% parameters fule feedback model
% Select feedback strength
% option 1= static
% option 2= weak-constant
% option 3= moderate-constant
% option 4= strong-constant
% option 5= weak-fading
% option 6= moderate-fading
% option 7= strong-fading
% Let's use option 6 or 7
option=1;

%area of forest extent from LANDFIRE ESP - this is the the km2 of forested
%land in the study area derived from Landfire environmental site potential
totalarea=77440;

% reference fuel memory (years)
nmem=30;

% annual max burned fraction ~ approximate 100-yr extreme value from the
% observed record
capannual=0.15;


%% load data
namefile = [dir_data,'fires/frap_forest_sierra_ncoast_year_1950_2021.mat'];
load(namefile,'FIRE')
%yor=(log(FIRE));
[C,IA,IB] = intersect(years_study,years_frap);
BA=FIRE(IB);
c=FIRE;

%climatic variables
filename = [dir_data,'nclimgrid/nclimgrid.mat'];
load(filename);
nclimgrid=nclimgrid(:,(11*12+1):end); %1971-2021
TSMAX_obs = zeros(length(years_study),1)*NaN;
for iyear=1:length(years_study)
    i1 = (iyear - 1) * 12 + best_start_tx;
    i2 = (iyear - 1) * 12 + best_stop_tx;
    TSMAX_obs(iyear) = mean(nclimgrid(1,i1:i2));
    PR_obs(iyear) = sum(nclimgrid(2,i1:i2));
end

%std over the baseline period 1995-2014
% [~,~,iok]=intersect(base_period, years_study);
% avg=nanmean(TSMAX_obs(iok),1);
% devstd=std(TSMAX_obs(iok),'omitnan');
% Z=(TSMAX_obs-avg)/devstd;

%% GCM

listScen ={'ssp245';'ssp585'};

namefile = [dir_data 'gcms/list_gcms_future.mat'];
load(namefile) %model_name

namefile = [dir_data 'gcms/tasmax-historical-1950-2014-monthly.mat'];
load(namefile) %historical
iok=find(~isnan(historical(1,:)));

namefile = [dir_data 'gcms/tasmax-ssp245-2015-2100-monthly.mat'];
load(namefile) %ssp245

namefile = [dir_data 'gcms/tasmax-ssp585-2015-2100-monthly.mat'];
load(namefile) %ssp585

m = matfile([dir_data 'gcms/pr-historical-1950-2014-monthly.mat']);
historical_p=m.historical;

m = matfile([dir_data 'gcms/pr-ssp245-2015-2100-monthly.mat']);
ssp245_p=m.ssp245;

m = matfile([dir_data 'gcms/pr-ssp585-2015-2100-monthly.mat']);
ssp585_p=m.ssp585;


model_name=model_name(iok);
historical=historical(:,iok);
ssp245=ssp245(:,iok);
ssp585=ssp585(:,iok);
historical_p=historical_p(:,iok);
ssp245_p=ssp245_p(:,iok);
ssp585_p=ssp585_p(:,iok);

model=cell(size(ssp585,2),1);
member=cell(size(ssp585,2),1);

for im=1:size(ssp585,2)
    aux=char(model_name(im));
    nchr=find(char(model_name(im)) == '_');
    model{im}=(aux(1:nchr-1));
    member{im}=aux(nchr+1:end);
end

models=unique(model);
length(models)

for imodel=1:length(models)
    models(imodel)
    imembers=find(strcmp(model,models(imodel)));
    members=member(imembers)
end



%figure;plot(ssp585)

%% ssp245 & ssp585
%[C,IA,IB] = intersect(years_sim_final,years_sim);
%[D,IC,ID] = intersect(years_sim_final,years_study);

%BCC-CSM2-MR_ssp245_r1i1p1f1
TSMAX_ssp245 = zeros([length(years_sim) size(models,1)] )*NaN;            % preallocation
TSMAX_ssp585 = TSMAX_ssp245; % preallocation
PR_ssp245 = zeros([length(years_sim) size(models,1)] )*NaN;            % preallocation
PR_ssp585 = TSMAX_ssp245; % preallocation
kk=0;
for imodel=1:length(models)
    models(imodel)
    imembers=find(strcmp(model,models(imodel)));
    members=member(imembers);
    for imember=1:length(members)
        kk=kk+1;
        members(imember)
        aux_ssp245=[historical(:,kk) ;ssp245(:,kk)];
        aux_ssp585=[historical(:,kk) ;ssp585(:,kk)];
        aux_ssp245_p=[historical_p(:,kk) ;ssp245_p(:,kk)];
        aux_ssp585_p=[historical_p(:,kk) ;ssp585_p(:,kk)];
        for iyear=1:length(years_sim)
            i1 = (iyear - 1) * 12 + best_start_tx;
            i2 = (iyear - 1) * 12 + best_stop_tx;
            TSMAX_245(iyear) = mean(aux_ssp245(i1:i2));
            TSMAX_585(iyear) = mean(aux_ssp585(i1:i2));
            PR_245(iyear) = mean(aux_ssp245_p(i1:i2));
            PR_585(iyear) = mean(aux_ssp585_p(i1:i2));
        end
        TSMAX_ssp245(:,kk)=TSMAX_245;
        TSMAX_ssp585(:,kk)=TSMAX_585;
        PR_ssp245(:,kk)=PR_245;
        PR_ssp585(:,kk)=PR_585;

    end
end
%figure;plot(cmip6_ssp245)
%figure;plot(cmip6_ssp585)

% bias correct output
nmod=size(TSMAX_ssp245,2);
[cmip6_ssp245]=biascorrect(TSMAX_obs,TSMAX_ssp245,nmod);
[cmip6_ssp585]=biascorrect(TSMAX_obs,TSMAX_ssp585,nmod);

[cmip6_ssp245_p]=biascorrect(PR_obs,PR_ssp245,nmod);
[cmip6_ssp585_p]=biascorrect(PR_obs,PR_ssp585,nmod);

for i=1:118
    cmip6_ssp585_t=modify_tp(cmip6_ssp585,cmip6_ssp585_p,1:71);
    cmip6_ssp245_t=modify_tp(cmip6_ssp245,cmip6_ssp245_p,1:71);
end

%%TSMAX changes
%simply averaging all ensemble members as ensemble mean would bias the result towards the models with more members. 
%That's why we average all members for each models before averaging all models

ino=find(strcmp(model,'CanESM5-CanOE') | strcmp(model,'CIESM') | strcmp(model,'GISS-E2-1-G'));
model_ok=model;
model_ok(ino)=[];
member_ok=member;
member_ok(ino)=[];
single_model_ok=unique(model_ok);
num_sim_tot=size(member_ok,1);

TSMAX_values_24_ssp245=zeros(size(cmip6_ssp245,1),length(single_model_ok))*NaN;
TSMAX_values_24_ssp585=zeros(size(cmip6_ssp245,1),length(single_model_ok))*NaN;
PR_values_24_ssp245=zeros(size(cmip6_ssp245,1),length(single_model_ok))*NaN;
PR_values_24_ssp585=zeros(size(cmip6_ssp245,1),length(single_model_ok))*NaN;
for im=1:length(single_model_ok)
    num_member=find(strcmp(model_ok,single_model_ok(im)));
    TSMAX_values_24_ssp245(:,im)=nanmean(cmip6_ssp245_t(:,num_member),2);
    TSMAX_values_24_ssp585(:,im)=nanmean(cmip6_ssp585_t(:,num_member),2);
    PR_values_24_ssp245(:,im)=nanmean(cmip6_ssp245_p(:,num_member),2);
    PR_values_24_ssp585(:,im)=nanmean(cmip6_ssp585_p(:,num_member),2);
end

% filename = [dir_data 'gcms/TSMAX_values_24_ssp245.mat'];
% save(filename,'TSMAX_values_24_ssp245','TSMAX_values_24_ssp245') 
% 
% filename = [dir_data 'gcms/TSMAX_values_24_ssp245.mat'];
% save(filename,'TSMAX_values_24_ssp585','TSMAX_values_24_ssp585') 

[~,~,Ipres] = intersect(base_period,years_sim);
[~,~,Ifut] = intersect(2031:2050,years_sim);
TSMAX_changes_ssp245(:)=mean(TSMAX_values_24_ssp245(Ifut,:),1)-mean(TSMAX_values_24_ssp245(Ipres,:),1);
prctile(TSMAX_changes_ssp245,[2.5,25,50,75,97.5])
TSMAX_changes_ssp585(:)=mean(TSMAX_values_24_ssp585(Ifut,:),1)-mean(TSMAX_values_24_ssp585(Ipres,:),1);
prctile(TSMAX_changes_ssp585,[2.5,25,50,75,97.5])

figure;plot(PR_values_24_ssp585)


PR_changes_ssp245(:)=100*(mean(PR_values_24_ssp245(Ifut,:),1)-mean(PR_values_24_ssp245(Ipres,:),1))./mean(PR_values_24_ssp245(Ipres,:),1);
prctile(PR_changes_ssp245,[2.5,25,50,75,97.5])
PR_changes_ssp585(:)=100*(mean(PR_values_24_ssp585(Ifut,:),1)-mean(PR_values_24_ssp585(Ipres,:),1))./mean(PR_values_24_ssp585(Ipres,:),1);
prctile(PR_changes_ssp585,[2.5,25,50,75,97.5])

%% fuel feedback model
wus = c(1:21,1); % this is data from 1950-1970
% I will recycle this three times to emulate fire suppression era data
wus=repmat(wus,[1 3]);
wus=wus(:);

TSMAX_obst=modify_tp(TSMAX_obs',PR_obs,1:51);

Z=(TSMAX_obst-mean(TSMAX_obst))./std(TSMAX_obst); %standardized over the 1971-2021 period
Z2=(PR_obs-mean(PR_obs))./std(PR_obs);

%standardized over the 1971-2021 period
cmip6_ssp245_t=(cmip6_ssp245_t-repmat(nanmean(cmip6_ssp245_t(22:71,:),1),[151 1]))./repmat(nanstd(cmip6_ssp245_t(22:71,:),[],1),[151 1]);
cmip6_ssp585_t=(cmip6_ssp585_t-repmat(nanmean(cmip6_ssp585_t(22:71,:),1),[151 1]))./repmat(nanstd(cmip6_ssp585_t(22:71,:),[],1),[151 1]);
cmip6_ssp245_p=(cmip6_ssp245_p-repmat(nanmean(cmip6_ssp245_p(22:71,:),1),[151 1]))./repmat(nanstd(cmip6_ssp245_p(22:71,:),[],1),[151 1]);
cmip6_ssp585_p=(cmip6_ssp585_p-repmat(nanmean(cmip6_ssp585_p(22:71,:),1),[151 1]))./repmat(nanstd(cmip6_ssp585_p(22:71,:),[],1),[151 1]);

%years_sim(1:71)
%years_sim(22:72)

% loop over cmip models
for imodel=1:nmod
    imodel

    % spin up historical fuel estimates based on fire activity
    [firehistory,rat1(:)]=spinupfirehistory(wus,option-1,totalarea,nmem);
    fire1950=firehistory(:,end-20);
    firehistory=firehistory(:,end);
    [FFA,ra,addr,raold,permold]=ba_feedbackmodel2(BA,Z,Z2,totalarea,cmip6_ssp245_t(:,imodel),cmip6_ssp245_p(:,imodel),option-1,capannual,nmem,firehistory,fire1950);
    % if model==1
    %    rat1(69:105)=raold+permold;
    % end
    ratg_ssp245(:,:,imodel)=ra;
    %rval(:,:,model)=r;
    burned_ssp245(:,:,imodel)=FFA;
    %arat(:,:,model)=addr;
    %oldrat(:,:)=raold;

    [FFA,ra,addr,raold,permold]=ba_feedbackmodel2(BA,Z,Z2,totalarea,cmip6_ssp585_t(:,imodel),cmip6_ssp585_p(:,imodel),option-1,capannual,nmem,firehistory,fire1950);
    % if model==1
    %    rat1(69:105)=raold+permold;
    % end
    ratg_ssp585(:,:,imodel)=ra;
    %rval(:,:,model)=r;
    burned_ssp585(:,:,imodel)=FFA;
    %arat(:,:,model)=addr;
    %oldrat(:,:)=raold;

end
% figure; plot(burned_ssp585(:,:,imodel))

%% calculate ensemble means

%simply averaging all ensemble members as ensemble mean would bias the result towards the models with more members. 
%That's why we average all members for each models before averaging all models

ino=find(strcmp(model,'CanESM5-CanOE') | strcmp(model,'CIESM') | strcmp(model,'GISS-E2-1-G'));
model_ok=model;
model_ok(ino)=[];
member_ok=member;
member_ok(ino)=[];
single_model_ok=unique(model_ok);
num_sim_tot=size(member_ok,1);


BA_values_24_ssp245=zeros(size(burned_ssp245,1),NB,length(single_model_ok))*NaN;
BA_values_24_ssp585=zeros(size(burned_ssp585,1),NB,length(single_model_ok))*NaN;
for im=1:length(single_model_ok)
    num_member=find(strcmp(model_ok,single_model_ok(im)));
    BA_values_24_ssp245(:,:,im)=nanmean(burned_ssp245(:,:,num_member),3);
    BA_values_24_ssp585(:,:,im)=nanmean(burned_ssp585(:,:,num_member),3);
end

filename = [dir_data 'gcms/BA_option_',num2str(option),'_tp.mat'];
save(filename,'BA_values_24_ssp245','BA_values_24_ssp585')  

%BA_values_24_ssp245(1:72,:,:)=repmat(BA_values_24_ssp245(1:72,1,1),[1 1000 length(single_model_ok)]);
%BA_values_24_ssp585(1:72,:,:)=repmat(BA_values_24_ssp585(1:72,1,1),[1 1000 length(single_model_ok)]);

BA_values_24_ssp245_reshaped=reshape(BA_values_24_ssp245,[size(BA_values_24_ssp245,1),size(BA_values_24_ssp245,2)*size(BA_values_24_ssp245,3)]);
BA_values_24_ssp585_reshaped=reshape(BA_values_24_ssp585,[size(BA_values_24_ssp585,1),size(BA_values_24_ssp585,2)*size(BA_values_24_ssp585,3)]);



figure; hold on;
%BA_values_24_ssp245_reshaped(1:72,:)=repmat(BA_values_24_ssp245_reshaped(1:72,1),[1  size(BA_values_24_ssp245_reshaped,2)]);
plot(1971:2100,movmedian(prctile(BA_values_24_ssp245_reshaped(22:151,:),[50],2),11),'k')
plot(1971:2100,movmedian(prctile(BA_values_24_ssp245_reshaped(22:151,:),[5],2),11),'--k')
plot(1971:2100,movmedian(prctile(BA_values_24_ssp245_reshaped(22:151,:),[95],2),11),'--k')

BA_values_24_ssp585_reshaped(1:72,:)=repmat(BA_values_24_ssp585_reshaped(1:72,1),[1  size(BA_values_24_ssp245_reshaped,2)]);
plot(1971:2100,movmedian(prctile(BA_values_24_ssp585_reshaped(22:151,:),[50],2),11),'r')
plot(1971:2100,movmedian(prctile(BA_values_24_ssp585_reshaped(22:151,:),[5],2),11),'--r')
plot(1971:2100,movmedian(prctile(BA_values_24_ssp585_reshaped(22:151,:),[95],2),11),'--r')
ylabel('Annual forest BA (km2) - 11-yr moving mean')

[~,~,ia5] = intersect(2031:2050,years_sim);
BA_values_24_ssp245_mod_cost=reshape(BA_values_24_ssp245,[size(BA_values_24_ssp245,1),size(BA_values_24_ssp245,2)*size(BA_values_24_ssp245,3)]);
100*(prctile(mean(BA_values_24_ssp245_mod_cost(ia5,:),1),[25 50 75])-2037)/2037