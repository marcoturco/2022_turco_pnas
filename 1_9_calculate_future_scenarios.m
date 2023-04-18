%% plot figure 2b, 2c, 2d
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
%addpath '~/Dropbox/estcena/scripts/fires_california/scripts/fuelfeedbackmodel_MT/'
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
end
%Z=(TSMAX_obs-avg)/devstd;
Z=(TSMAX_obs-mean(TSMAX_obs))./std(TSMAX_obs); %standardized over the 1971-2021 period

%% GCM

listScen ={'ssp245';'ssp585'};

namefile = [dir_data 'gcms/list_gcms_future.mat'];
load(namefile) %model_name

namefile = [dir_data 'gcms/tasmax-historical-1950-2014-monthly.mat'];
load(namefile) %historical
iok=find(~isnan(historical(1,:)));

namefile = [dir_data 'gcms/tasmax-ssp245-2015-2100-monthly.mat'];
load(namefile) %ssp245

namefile = [dir_data,'gcms/tasmax-ssp585-2015-2100-monthly.mat'];
load(namefile) %ssp585

model_name=model_name(iok);
historical=historical(:,iok);
ssp245=ssp245(:,iok);
ssp585=ssp585(:,iok);

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
        for iyear=1:length(years_sim)
            i1 = (iyear - 1) * 12 + best_start_tx;
            i2 = (iyear - 1) * 12 + best_stop_tx;
            TSMAX_245(iyear) = mean(aux_ssp245(i1:i2));
            TSMAX_585(iyear) = mean(aux_ssp585(i1:i2));
        end
        TSMAX_ssp245(:,kk)=TSMAX_245;
        TSMAX_ssp585(:,kk)=TSMAX_585;
        %TSMAX_245=TSMAX_245(IB);
        %TSMAX_585=TSMAX_585(IB);
        %
        %         % bias adjustment
        %         [aux1,aux2]=bias_correction_ecdf...
        %             (TSMAX_obs,TSMAX_245(ID),TSMAX_245);
        %         TSMAX_anom=scale_base_period(aux2',base_period,years_study);
        %         TSMAX_ssp245(:,kk)=TSMAX_anom;
        %         % bias adjustment
        %         [aux1,aux2]=bias_correction_ecdf...
        %             (TSMAX_obs,TSMAX_585(ID),TSMAX_585);
        %         TSMAX_anom=scale_base_period(aux2',base_period,years_study);
        %         TSMAX_ssp585(:,kk)=TSMAX_anom;
        %         %figure;plot(aux2)
    end
end
%figure;plot(cmip6_ssp245)
%figure;plot(cmip6_ssp585)

% bias correct output
nmod=size(TSMAX_ssp245,2);
[cmip6_ssp245]=biascorrect(TSMAX_obs,TSMAX_ssp245,nmod);
[cmip6_ssp585]=biascorrect(TSMAX_obs,TSMAX_ssp585,nmod);

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
for im=1:length(single_model_ok)
    num_member=find(strcmp(model_ok,single_model_ok(im)));
    TSMAX_values_24_ssp245(:,im)=nanmean(cmip6_ssp245(:,num_member),2);
    TSMAX_values_24_ssp585(:,im)=nanmean(cmip6_ssp585(:,num_member),2);
end

filename = [dir_data 'gcms/TSMAX_values_24_ssp245.mat'];
save(filename,'TSMAX_values_24_ssp245','TSMAX_values_24_ssp245') 

filename = [dir_data 'gcms/TSMAX_values_24_ssp245.mat'];
save(filename,'TSMAX_values_24_ssp585','TSMAX_values_24_ssp585') 

[~,~,Ipres] = intersect(base_period,years_sim);
[~,~,Ifut] = intersect(2031:2050,years_sim);
TSMAX_changes_ssp245(:)=mean(TSMAX_values_24_ssp245(Ifut,:),1)-mean(TSMAX_values_24_ssp245(Ipres,:),1);
prctile(TSMAX_changes_ssp245,[2.5,25,50,75,97.5])
TSMAX_changes_ssp585(:)=mean(TSMAX_values_24_ssp585(Ifut,:),1)-mean(TSMAX_values_24_ssp585(Ipres,:),1);
prctile(TSMAX_changes_ssp585,[2.5,25,50,75,97.5])

%% fuel feedback model
wus = c(1:21,1); % this is data from 1950-1970
% I will recycle this three times to emulate fire suppression era data
wus=repmat(wus,[1 3]);
wus=wus(:);




model_name(imodel)
% loop over cmip models
for imodel=1:nmod
    imodel
    % spin up historical fuel estimates based on fire activity
    [firehistory,rat1(:)]=spinupfirehistory(wus,option-1,totalarea,nmem);
    fire1950=firehistory(:,end-20);
    firehistory=firehistory(:,end);

    avg_tsmax=nanmean(cmip6_ssp245(22:71,imodel),1);
    devstd_tsmax=std(cmip6_ssp245(22:71,imodel),'omitnan');
    futureZ=(cmip6_ssp245(:,imodel)-avg_tsmax)/devstd_tsmax;
    [r,p,FFA,ra,addr,raold,permold]=ba_feedbackmodel(BA,Z,totalarea,futureZ,option-1,capannual,nmem,firehistory,fire1950);
    % if model==1
    %    rat1(69:105)=raold+permold;
    % end
    ratg_ssp245(:,:,imodel)=ra;
    %rval(:,:,model)=r;
    burned_ssp245(:,:,imodel)=FFA;
    %arat(:,:,model)=addr;
    %oldrat(:,:)=raold;

    % spin up historical fuel estimates based on fire activity
    [firehistory,rat1(:)]=spinupfirehistory(wus,option-1,totalarea,nmem);
    fire1950=firehistory(:,end-20);
    firehistory=firehistory(:,end); % if model==1
    
    avg_tsmax=nanmean(cmip6_ssp585(22:71,imodel),1);
    devstd_tsmax=std(cmip6_ssp585(22:71,imodel),'omitnan');
    futureZ=(cmip6_ssp585(:,imodel)-avg_tsmax)/devstd_tsmax;
    [r,p,FFA,ra,addr,raold,permold]=ba_feedbackmodel(BA,Z,totalarea,futureZ,option-1,capannual,nmem,firehistory,fire1950);
    
    %    rat1(69:105)=raold+permold;
    % end
    ratg_ssp585(:,:,imodel)=ra;
    %rval(:,:,model)=r;
    burned_ssp585(:,:,imodel)=FFA;
    %arat(:,:,model)=addr;
    %oldrat(:,:)=raold;

end

figure; plot(ra), ylim([0 1])

figure; plot(burned_ssp585(:,:,imodel))

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

filename = [dir_data 'gcms/BA_option_',num2str(option),'.mat'];
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
