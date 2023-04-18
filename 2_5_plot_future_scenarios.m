%% plot figure 3
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
version='';

%%
years_study=1971:2021;
base_period=1995:2014;
years_frap=1950:2021;
years_hist = 1950:2014; %years hist-nat simulations
years_fut = 2015:2100; %years future simulations
years_sim=1950:2100;
years_sim_final=1971:2050;
col_mod_obs_line=[110 110 110]/255;
[~,~,ia1] = intersect(1971:2050,years_sim);
[~,~,ia2] = intersect(1971:2014,years_sim);
[~,~,ia3] = intersect(2014:2050,years_sim);
[~,~,ia4] = intersect(1971:2050,years_sim);

[~,~,ia6] = intersect(base_period,years_study);
font_size=12;

%% frap 
namefile = [dir_data,'fires/frap_forest_sierra_ncoast_year.mat'];
load(namefile,'FIRE')

filename = [dir_data 'gcms/BA_option_1',version,'.mat']; %static model
load(filename) %,'BA_values_24_ssp245','BA_values_24_ssp585') 
BA=FIRE;

BA_values_24_ssp245_static=reshape(BA_values_24_ssp245,[size(BA_values_24_ssp245,1),size(BA_values_24_ssp245,2)*size(BA_values_24_ssp245,3)]);
BA_values_24_ssp585_static=reshape(BA_values_24_ssp585,[size(BA_values_24_ssp585,1),size(BA_values_24_ssp585,2)*size(BA_values_24_ssp585,3)]);

size(BA_values_24_ssp585_static)


% [~,~,Ipres] = intersect(base_period,years_sim);
% [~,~,Ifut] = intersect(2031:2050,years_sim);
% size(BA)
% for i=1:size(BA_values_24_ssp585_static,2)
%     BA_changes_ssp245(i)=(mean(BA_values_24_ssp245_static(Ifut,i),1)-mean(BA_values_24_ssp245_static(Ipres,i),1))./mean(BA_values_24_ssp245_static(Ipres,i),1);
%     BA_changes_ssp585(i)=(mean(BA_values_24_ssp585_static(Ifut,i),1)-mean(BA_values_24_ssp585_static(Ipres,i),1))./mean(BA_values_24_ssp585_static(Ipres,i),1);
% end
% prctile(BA_changes_ssp245,[25,50,75])
% prctile(BA_changes_ssp585,[25,50,75])


aux2(:,1)=movmean(median(BA_values_24_ssp585_static,2),21);
aux2(:,2)=movmean(prctile(BA_values_24_ssp585_static,[97.5],2),21);
aux2(:,3)=movmean(prctile(BA_values_24_ssp585_static,[2.5],2),21);

aux1(:,1)=movmean(median(BA_values_24_ssp245_static,2),21);
aux1(:,2)=movmean(prctile(BA_values_24_ssp245_static,[97.5],2),21);
aux1(:,3)=movmean(prctile(BA_values_24_ssp245_static,[2.5],2),21);

indices{1}=[1:size(BA_values_24_ssp245_static,2)];
boundary=cell(1,1);boundary{1}={'prc5','prc95'}; 

%% static
figure; hold on;
bar(years_study,BA,'r')
drawSpread(aux1(ia2,:),'xvalues',years_sim(ia2)','colorsg',[191 191 191]/255,'lines','no','alphasg',[0.5 0.5])
plot(years_sim(ia2),aux1(ia2,1),'color',[0 0 0]/255)
%drawSpread(aux1(ia3,:),'xvalues',years_sim(ia3)','colorsg',[255 165 0]/255,'lines','no','alphasg',[0.5 0.5])
%plot(years_sim(ia3),aux1(ia3,1),'color',[155 140 0]/255)
%drawSpread(aux2(ia3,:),'xvalues',years_sim(ia3)','colorsg',[243 209 213]/255,'lines','no','alphasg',[0.5 0.5])
%plot(years_sim(ia3),aux2(ia3,1),'color',[153 43 62]/255)
drawSpread(aux1(ia3,:),'xvalues',years_sim(ia3)','colorsg',[191 191 191]/255,'lines','no','alphasg',[0.5 0.5])
plot(years_sim(ia3),aux1(ia3,1),'color',[0 0 0]/255)
%legend('Obs.','95% CI CMIP6','Median CMIP6','Location','Best','AutoUpdate','off')
legend('Obs.','95% CI CMIP6','Median CMIP6','Location','northwest','AutoUpdate','off')
%legend boxoff
bar(years_study,BA,'r')
xlabel('Years','FontSize',font_size) 
ylabel('BA (km^2)','FontSize',font_size) 
%gridxy(years(ith),'Color','k','Linestyle',':') ;
set(gca,'FontSize',font_size,'xlim',[1970 2051])
file=[dir_out,'fig3a.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),version,'.pdf'] ,'-painters')

%%  feedbacks plot
filename = [dir_data 'gcms/BA_option_3',version,'.mat']; 
% option 3= moderate-constant
% option 6= moderate-fading
load(filename) %,'BA_values_24_ssp245','BA_values_24_ssp585') 
% BA_values_24_ssp245(1:72,:,:)=repmat(BA_values_24_ssp245(1:72,1,1),[1 1000 24]);
% BA_values_24_ssp585(1:72,:,:)=repmat(BA_values_24_ssp585(1:72,1,1),[1 1000 24]);
BA_values_24_ssp245_mod_cost=reshape(BA_values_24_ssp245,[size(BA_values_24_ssp245,1),size(BA_values_24_ssp245,2)*size(BA_values_24_ssp245,3)]);
BA_values_24_ssp585_mod_cost=reshape(BA_values_24_ssp585,[size(BA_values_24_ssp585,1),size(BA_values_24_ssp585,2)*size(BA_values_24_ssp585,3)]);
BA_values_24_ssp245_mod_cost(1:72,:)=BA_values_24_ssp245_static(1:72,:);
BA_values_24_ssp585_mod_cost(1:72,:)=BA_values_24_ssp585_static(1:72,:);

filename = [dir_data 'gcms/BA_option_6',version,'.mat']; 
% option 3= moderate-constant
% option 6= moderate-fading
load(filename) %,'BA_values_24_ssp245','BA_values_24_ssp585') 
% BA_values_24_ssp245(1:72,:,:)=repmat(BA_values_24_ssp245(1:72,1,1),[1 1000 24]);
% BA_values_24_ssp585(1:72,:,:)=repmat(BA_values_24_ssp585(1:72,1,1),[1 1000 24]);
BA_values_24_ssp245_mod_fad=reshape(BA_values_24_ssp245,[size(BA_values_24_ssp245,1),size(BA_values_24_ssp245,2)*size(BA_values_24_ssp245,3)]);
BA_values_24_ssp585_mod_fad=reshape(BA_values_24_ssp585,[size(BA_values_24_ssp585,1),size(BA_values_24_ssp585,2)*size(BA_values_24_ssp585,3)]);
BA_values_24_ssp245_mod_fad(1:72,:)=BA_values_24_ssp245_static(1:72,:);
BA_values_24_ssp585_mod_fad(1:72,:)=BA_values_24_ssp585_static(1:72,:);


aux3(:,1)=movmean(median(BA_values_24_ssp245_mod_cost,2),21);
aux4(:,1)=movmean(median(BA_values_24_ssp585_mod_cost,2),21);
aux5(:,1)=movmean(median(BA_values_24_ssp245_mod_fad,2),21);
aux6(:,1)=movmean(median(BA_values_24_ssp585_mod_fad,2),21);

figure;
plot(years_study,movmean(BA,21),'k')
hold on
plot(years_sim(ia4),aux1(ia4,1),'color',[191 191 191]/255)
plot(years_sim(ia4),aux2(ia4,1),'color',[150 150 150]/255)
plot(years_sim(ia4),aux3(ia4,1),'color',[210 105 30]/255)
plot(years_sim(ia4),aux4(ia4,1),'color',[165 42 42]/255)
plot(years_sim(ia4),aux5(ia4,1),'color',[170 189 216]/255)
plot(years_sim(ia4),aux6(ia4,1),'color',[78 143 188]/255)
legend('Observed','Static ssp245','Static ssp585','Moderate-Constant ssp245','Moderate-Constant ssp585','Moderate-Fading ssp245','Moderate-Fading ssp585','Location','Best')
legend boxoff
xlabel('Years','FontSize',font_size) 
ylabel('BA (km^2)','FontSize',font_size) 
%gridxy(years(ith),'Color','k','Linestyle',':') ;
set(gca,'FontSize',font_size,'xlim',[1970 2051])
file=[dir_out,'fig3b.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),version,'.pdf'] ,'-painters')


%% boxplot

dataplot=zeros(1,size(BA_values_24_ssp585_static,2),6)*NaN;
[~,~,ia5] = intersect(2031:2050,years_sim);
dataplot(1,1:size(BA_values_24_ssp245_static,2),1)=(mean(BA_values_24_ssp245_static(ia5,:),1));
dataplot(1,1:size(BA_values_24_ssp585_static,2),2)=(mean(BA_values_24_ssp585_static(ia5,:),1));
dataplot(1,1:size(BA_values_24_ssp245_mod_cost,2),3)=(mean(BA_values_24_ssp245_mod_cost(ia5,:),1));
dataplot(1,1:size(BA_values_24_ssp585_mod_cost,2),4)=(mean(BA_values_24_ssp585_mod_cost(ia5,:),1));
dataplot(1,1:size(BA_values_24_ssp245_mod_fad,2),5)=(mean(BA_values_24_ssp245_mod_fad(ia5,:),1));
dataplot(1,1:size(BA_values_24_ssp585_mod_fad,2),6)=(mean(BA_values_24_ssp585_mod_fad(ia5,:),1));
c=[191 191 191;150 150 150;210 105 30;165 42 42;170 189 216;78 143 188]/255;
% 

prctile(mean(BA_values_24_ssp245_static(ia5,:),1),[2.5 25 50 75 97.5])
prctile(mean(BA_values_24_ssp585_static(ia5,:),1),[2.5 25 50 75 97.5])

prctile(mean(BA_values_24_ssp245_mod_fad(ia5,:),1),[2.5 25 50 75 97.5])
prctile(mean(BA_values_24_ssp245_mod_cost(ia5,:),1),[2.5 25 50 75 97.5])

prctile(mean(BA_values_24_ssp585_mod_fad(ia5,:),1),[2.5 25 50 75 97.5])
prctile(mean(BA_values_24_ssp585_mod_cost(ia5,:),1),[2.5 25 50 75 97.5])


100*(prctile(mean(BA_values_24_ssp245_static(ia5,:),1),[25 50 75])-2037)/2037
100*(prctile(mean(BA_values_24_ssp585_static(ia5,:),1),[25 50 75])-2037)/2037
100*(prctile(mean(BA_values_24_ssp245_mod_fad(ia5,:),1),[25 50 75])-2037)/2037
100*(prctile(mean(BA_values_24_ssp585_mod_fad(ia5,:),1),[25 50 75])-2037)/2037
100*(prctile(mean(BA_values_24_ssp245_mod_cost(ia5,:),1),[25 50 75])-2037)/2037
100*(prctile(mean(BA_values_24_ssp585_mod_cost(ia5,:),1),[25 50 75])-2037)/2037



%mean(BA_values_24_ssp245_static(ia5,:),1);
figure; hold on;
aboxplot3(dataplot,'colormap',c,'colorrev',1)
ylim([mean(BA(ia6))*0.8 BA(end)*1.2])
gridxy([],mean(BA(ia6)))
text(0.55,mean(BA(ia6))+mean(BA(ia6))*0.3,['Obs. mean (1995-2014): ',num2str(mean(BA(ia6)),'%4.0f') ' km^2'])
gridxy([],mean(BA(32:51)))
text(0.55,mean(BA(32:51))-200,['Obs. mean (2002-2021): ',num2str(((mean(BA(32:51)))),'%4.0f') ' km^2'])
gridxy([],((BA(end))))
text(0.55,(round(BA(end)))+300,['Obs. 2021: ',num2str((((BA(end)))),'%4.0f') ' km^2'])
ylabel('Burned Area (km^2)','FontSize',font_size);
%legend('Observed','Static ssp245','Static ssp585','Moderate-Constant ssp245','Moderate-Constant ssp585','Moderate-Fading ssp245','Moderate-Fading ssp585','Location','Best')
ylim([700 10500])
xlabel_boxplot(1,:)=['Static ssp245           '];
xlabel_boxplot(2,:)=['Static ssp585           '];
xlabel_boxplot(3,:)=['Moderate-Constant ssp245'];
xlabel_boxplot(4,:)=['Moderate-Constant ssp585'];
xlabel_boxplot(5,:)=['Moderate-Fading ssp245  '];
xlabel_boxplot(6,:)=['Moderate-Fading ssp585  '];
set(gca,'XTick',1:size(xlabel_boxplot,1));
set(gca,'XTickLabel',xlabel_boxplot,'Fontsize',font_size);
set(gca,'FontSize',font_size)
set(gca,'FontSize',font_size)
nomeout=[dir_out,'fig3c',version,'.eps'];
print( gcf, '-depsc2', nomeout ,'-painters')
