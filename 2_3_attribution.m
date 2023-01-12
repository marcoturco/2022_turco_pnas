%% 
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years_study=1971:2021;
years_nat_all=1850:2021;
years_nat=1950:2021;
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
font_size=12;
base_period=1960:1982;
best_start_tx=4;
best_stop_tx=10;
alpha = 0.05;
NB=10000;
n=1; %linear trend

%% load obs climate and fires
%frap 
%frap 
namefile = [dir_data,'fires/frap_forest_sierra_ncoast_year.mat'];
load(namefile,'FIRE')
yor=(log(FIRE));


%climatic variables
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
[C,IA,IB] = intersect(years_study,years_climate);
TSMAX_obs_anom=TSMAX_obs_anom(IB);
TSMAX_obs=TSMAX_obs(IB);
mean(TSMAX_obs(31:50))-mean(TSMAX_obs(1:20))

%% load data hist-nat vs all
nomefile=[dir_data 'gcms/GCMs_attribution.mat'];
load(nomefile)
%figure;plot(TSMAX_ALL)
%figure;plot(TMAX_NAT)
%simply averaging all ensemble members as ensemble mean would bias the result towards the models with more members. 
%That?s why we average all members for each models before averaging all models



%iok2=sum(TSMAX_ALL,1)
%length(find(isnan(iok2)))
aux=sum(TSMAX_NAT,1);
iok=(find(~isnan(aux)));
aux=cell(1,length(iok));
for ii = 1:length(iok)
  aux{ii} = model_com{iok(ii)};
end
model_com=aux;
model_member_num=model_member_num(iok,:);
% 
% for im=1:model_member_num(end,1)
%     num_member=find(model_member_num(:,1)==im);
%     model_com{num_member}
%     length(num_member)
%     pause
% end

TSMAX_ALL=TSMAX_ALL(:,iok);
TSMAX_NAT=TSMAX_NAT(:,iok);
TSMAX_ALL_ANOM=TSMAX_ALL_ANOM(:,iok);
TSMAX_NAT_ANOM=TSMAX_NAT_ANOM(:,iok);

for im=1:model_member_num(end,1)
    num_member=find(model_member_num(:,1)==im);
    TSMAX_ALL_12(:,im)=nansum(TSMAX_ALL(:,num_member),2)/length(num_member);
    TSMAX_NAT_12(:,im)=nansum(TSMAX_NAT(:,num_member),2)/length(num_member);
    TSMAX_ALL_12_ANOM(:,im)=nansum(TSMAX_ALL_ANOM(:,num_member),2)/length(num_member);
    TSMAX_NAT_12_ANOM(:,im)=nansum(TSMAX_NAT_ANOM(:,num_member),2)/length(num_member);
end


for im=1:model_member_num(end,1)
    num_member=find(model_member_num(:,1)==im);
    TSMAX_ALL_1(:,im)=TSMAX_ALL(:,num_member(1));
    TSMAX_NAT_1(:,im)=TSMAX_NAT(:,num_member(1));
    TSMAX_ALL_1_ANOM(:,im)=TSMAX_ALL_ANOM(:,num_member(1));
    TSMAX_NAT_1_ANOM(:,im)=TSMAX_NAT_ANOM(:,num_member(1));
end

figure;plot(TSMAX_NAT_1)
figure;plot(TSMAX_NAT_1_ANOM)


figure;plot(TSMAX_NAT_12)
figure;plot(TSMAX_NAT_12_ANOM)

size(TSMAX_NAT_1_ANOM)
TSMAX_NAT_1_ANOM=TSMAX_NAT_1_ANOM(IB,:);
TSMAX_ALL_1_ANOM=TSMAX_ALL_1_ANOM(IB,:);
TSMAX_NAT_12_ANOM=TSMAX_NAT_12_ANOM(IB,:);
TSMAX_ALL_12_ANOM=TSMAX_ALL_12_ANOM(IB,:);
TSMAX_NAT_ANOM=TSMAX_NAT_ANOM(IB,:);
TSMAX_ALL_ANOM=TSMAX_ALL_ANOM(IB,:);

%% regression model with observations
Xorig1 = [ones(size(TSMAX_obs,1),1)  TSMAX_obs_anom ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
bootdat1 = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
%bootb1 = sort(bootb1);
%bootCI1 = [bootb1(25,:); bootb1(975,:)]  %Report 5% confidence intervals
B1

sigma=std(R);
oos1=[]; 
for iboot=1:NB
    noise=normrnd(0,sigma,length(yor),1);
    Y = (Xorig1*B1)+noise;
    oos1=[oos1,Y];
end

for i=1:length(years_study)
    Ymed(i)=prctile(exp(oos1(i,:)),50);
    b1(i,:)=[prctile(exp(oos1(i,:)),2.5),prctile(exp(oos1(i,:)),97.5)];
end 

%% TSMAX AVERAGED MEMBERS FOR EACH MODEL

trend_obs=movmean(TSMAX_obs_anom,21);
trend_nat=movmean(nanmean(TSMAX_NAT_12_ANOM,2),21); 
trend_all=movmean(nanmean(TSMAX_ALL_12_ANOM,2),21); 

indices{1}=[1:size(TSMAX_ALL_12_ANOM,2)];
boundary=cell(1,1);boundary{1}={'prc2.5','prc97.5'}; 
col_all_spread=[237 220 202]/255;
col_all_line=[130 90 38]/255;
col_nat_spread=[181 209 206]/255;
col_nat_line=[65 116 110]/255;

figure; hold on
%plot(years_study,prctile(TSMAX_ALL,[50],2),'-','LineWidth',2,'Color',col_all_line)
%plot(years_study,prctile(TSMAX_NAT,[50],2),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,trend_all,'-','LineWidth',2,'Color',col_all_line)
plot(years_study,trend_nat,'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,trend_obs,'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
%plot(years_study,running_avg(TSMAX_obs,num_years_runmean),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
legend('CMIP6-ALL','CMIP6-NAT','Observed','Location','Best','AutoUpdate','off')
legend boxoff
%drawSpread(TSMAX_NAT,'xvalues',years_study','indexesg',indices,'colorsg',col_nat_spread,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]);
%drawSpread(TSMAX_ALL,'xvalues',years_study','indexesg',indices,'colorsg',col_all_spread,'lines','no','boundary',boundary,'alphasg',[0.5 0.5])
plot(years_study,TSMAX_NAT_12_ANOM,'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,TSMAX_ALL_12_ANOM,'-','LineWidth',0.5,'Color',col_all_spread)
plot(years_study,(TSMAX_obs_anom(:)'),'k:o','LineWidth',1)
plot(years_study,trend_all,'-','LineWidth',2,'Color',col_all_line)
plot(years_study,trend_nat,'-','LineWidth',2,'Color',col_nat_line)
%plot(years_study,running_avg(TSMAX_obs,num_years_runmean),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(years_study,(TSMAX_obs_anom(:)'),'k:o','LineWidth',1)
plot(years_study,trend_obs,'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('years','FontSize',font_size) 
ylabel('TS_{max} [C^\circ]','FontSize',font_size) 
set(gca,'FontSize',font_size,'xlim',[years_study(1)-1 years_study(end)+1])
%gridxy([],[-2.5:0.5:2.5],'Color',[121 121 121]/255,'Linestyle','-');
%ylim([-2.5 2.5])
file=[dir_out,'figure_2a_tsmax_all_members.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')




mean(TSMAX_obs_anom(26:51))-mean(TSMAX_obs_anom(1:25))

prctile((mean(TSMAX_NAT_12_ANOM(26:51,:),1)-mean(TSMAX_NAT_12_ANOM(1:25,:),1)),[2.5 50 97.5],2)
prctile((mean(TSMAX_ALL_12_ANOM(26:51,:),1)-mean(TSMAX_ALL_12_ANOM(1:25,:),1)),[2.5 50 97.5],2)

prctile((mean(TSMAX_NAT_1_ANOM(26:51,:),1)-mean(TSMAX_NAT_1_ANOM(1:25,:),1)),[2.5 50 97.5],2)
prctile((mean(TSMAX_ALL_1_ANOM(26:51,:),1)-mean(TSMAX_ALL_1_ANOM(1:25,:),1)),[2.5 50 97.5],2)


%% BA
k=0;
for ifile=1:model_member_num(end,1)
    for ib=1:NB
        k=k+1;
        noise=normrnd(0,sigma,length(yor),1);
        pred_all_1(:,ifile,ib) = B1(1)  + B1(2) * TSMAX_ALL_1_ANOM(:,ifile) + noise;
        pred_nat_1(:,ifile,ib) = B1(1)  + B1(2) * TSMAX_NAT_1_ANOM(:,ifile) + noise;
    end
end

pred_all_1_plot(:,:)=nanmean(pred_all_1,3);
pred_nat_1_plot(:,:)=nanmean(pred_nat_1,3);


k=0;
for ifile=1:size(TSMAX_ALL_ANOM,2)
    for ib=1:NB
        noise=normrnd(0,sigma,length(yor),1);
        pred_all_all(:,ifile,ib) = B1(1)  + B1(2) * TSMAX_ALL_ANOM(:,ifile) + noise;
        pred_nat_all(:,ifile,ib) = B1(1)  + B1(2) * TSMAX_NAT_ANOM(:,ifile) + noise;        
    end
end

for im=1:model_member_num(end,1)
    num_member=find(model_member_num(:,1)==im);
    %i1=(num_member(1)-1)*NB+1;
    %i2=(num_member(end))*NB;
    pred_all_12(:,im,:)=nanmean(pred_all_all(:,num_member,:),2);
    pred_nat_12(:,im,:)=nanmean(pred_nat_all(:,num_member,:),2);
end

pred_all_12_plot(:,:)=nanmean(pred_all_12,3);
pred_nat_12_plot(:,:)=nanmean(pred_nat_12,3);

%% BA PLOT AVERAGED MEMBERS
%plot with ensemble mean
%trend_obs=polyval(polyfit(years_study,(yor(:))',1),years_study);
%trend_nat=polyval(polyfit(years_study,(prctile((pred_nat_12_plot),[50],2))',1),years_study);
%trend_all=polyval(polyfit(years_study,(prctile((pred_all_12_plot),[50],2))',1),years_study);
%trend_obs=(yor);
%trend_nat=(nanmean(pred_nat_12_plot,2)); 
%trend_all=(nanmean(pred_all_12_plot,2)); 
trend_obs=movmean(yor,21);
trend_nat=movmean(nanmean(pred_nat_12_plot,2),21); 
trend_all=movmean(nanmean(pred_all_12_plot,2),21); 


% Ba_orig=exp(yor)/1000;
% [sigBA,BAe]=testTrend_MT_sen(Ba_orig);
% TBAe=100*(BAe.*(length(Ba_orig)))./nanmean(Ba_orig,1)
% BAe.*(length(Ba_orig))*1000
% mean(Ba_orig(1:29))
% mean(Ba_orig(30:50))
% mean(Ba_orig(30:50))/mean(Ba_orig(1:29))
% 
% mean(Ba_orig(1:25))
% mean(Ba_orig(26:50))
% mean(Ba_orig(26:50))/mean(Ba_orig(1:25))

indices{1}=[1:size(pred_all_12_plot,2)];
boundary=cell(1,1);boundary{1}={'prc2.5','prc97.5'}; 

figure; hold on
plot(years_study,exp(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,exp(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,exp(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
legend('CMIP6-ALL','CMIP6-NAT','Observed','Location','Best','AutoUpdate','off')
legend boxoff
plot(years_study,exp(pred_nat_12_plot),'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,exp(pred_all_12_plot),'-','LineWidth',0.5,'Color',col_all_spread)
plot(years_study,exp(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,exp(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,exp(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(years_study,exp(yor(:)'),'k:o','LineWidth',1)
set(gca, 'YScale', 'log')
xlabel('years','FontSize',font_size) 
ylabel('Burned Area [km^2]','FontSize',font_size) 
set(gca,'FontSize',font_size,'xlim',[years_study(1)-1 years_study(end)+1])
%gridxy([],[-2.5:0.5:2.5],'Color',[121 121 121]/255,'Linestyle','-');
%ylim([-2.5 2.5])
file=[dir_out,'figure_2a_avg_members.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')

%%
trend_obs=movmean(yor,21);
trend_nat=movmean(nanmean(pred_nat_1_plot,2),21); 
trend_all=movmean(nanmean(pred_all_1_plot,2),21); 

indices{1}=[1:size(pred_all_1_plot,2)];
boundary=cell(1,1);boundary{1}={'prc2.5','prc97.5'}; 

figure; hold on
plot(years_study,exp(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,exp(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,exp(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
legend('CMIP6-ALL','CMIP6-NAT','Observed','Location','Best','AutoUpdate','off')
legend boxoff
plot(years_study,exp(pred_nat_1_plot),'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,exp(pred_all_1_plot),'-','LineWidth',0.5,'Color',col_all_spread)
plot(years_study,exp(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,exp(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,exp(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(years_study,exp(yor(:)'),'k:o','LineWidth',1)
set(gca, 'YScale', 'log')
xlabel('years','FontSize',font_size) 
ylabel('Burned Area [km^2]','FontSize',font_size) 
set(gca,'FontSize',font_size,'xlim',[years_study(1)-1 years_study(end)+1])
%gridxy([],[-2.5:0.5:2.5],'Color',[121 121 121]/255,'Linestyle','-');
%ylim([-2.5 2.5])
file=[dir_out,'figure_2a_ba_one_member.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')



%% regression D&A
%figure;plot(pred_nat_12_plot)
trend_nat=(nanmean(pred_nat_12_plot,2)); 
trend_all=(nanmean(pred_all_12_plot,2)); 


Xall = [ones(size(trend_all,1),1) trend_all]; %osservati
Xnat = [ones(size(trend_nat,1),1) trend_nat]; %osservati
for iy=1:22
    %[b,bint] = regress(yor(1:(29+iy)),Xall(1:(29+iy),:));
    %Ball(iy,1)=b(2);
    %Ball(iy,2:3)=bint(2,:);
    bootdat = [yor(1:(29+iy)),Xall(1:(29+iy),:)];
    aux= bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
    Ball_boot(iy,:,:)=aux;
    %[b,bint]=regress(yor(1:(29+iy)),Xnat(1:(29+iy),:));
    %Bnat(iy,1)=b(2);
    %Bnat(iy,2:3)=bint(2,:);
    bootdat = [yor(1:(29+iy)),Xnat(1:(29+iy),:)];
    aux=bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
    Bnat_boot(iy,:,:) = aux;
end

indices{1}=[1:size(Ball_boot,2)];
boundary=cell(1,1);boundary{1}={'prc5','prc95'}; 
col_all_spread=[237 220 202]/255;
col_all_line=[130 90 38]/255;
col_nat_spread=[181 209 206]/255;
col_nat_line=[65 116 110]/255;

figure; hold on
plot(2000:2021,prctile(Ball_boot(:,:,2),[50],2),'-','LineWidth',2,'Color',col_all_line)
plot(2000:2021,prctile(Bnat_boot(:,:,2),[50],2),'-','LineWidth',2,'Color',col_nat_line)
legend('CMIP6-ALL','CMIP6-NAT','Location','Best','AutoUpdate','off')
legend boxoff
drawSpread(Bnat_boot(:,:,2),'xvalues',(2000:2021)','indexesg',indices,'colorsg',col_nat_spread,'lines','no','boundary',boundary,'alphasg',[0.5 0.5])
drawSpread(Ball_boot(:,:,2),'xvalues',(2000:2021)','indexesg',indices,'colorsg',col_all_spread,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]);
xlabel('Last year of the serie','FontSize',18) 
ylabel('Scaling factor','FontSize',18) 
set(gca,'FontSize',18,'xlim',[2000,2021])
gridxy([],[0],'Color',[0 0 0]/255,'Linestyle','-','LineWidth',2);
ylim([-2.5 2.5])
file=[dir_out,'scaling_factor.eps']
% print( gcf, '-depsc2', file ,'-painters')
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')

%% regression D&A 1
%figure;plot(pred_nat_12_plot)
trend_nat=(nanmean(pred_nat_1_plot,2)); 
trend_all=(nanmean(pred_all_1_plot,2)); 


Xall = [ones(size(trend_all,1),1) trend_all]; %osservati
Xnat = [ones(size(trend_nat,1),1) trend_nat]; %osservati
for iy=1:22
    %[b,bint] = regress(yor(1:(29+iy)),Xall(1:(29+iy),:));
    %Ball(iy,1)=b(2);
    %Ball(iy,2:3)=bint(2,:);
    bootdat = [yor(1:(29+iy)),Xall(1:(29+iy),:)];
    aux= bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
    Ball_boot(iy,:,:)=aux;
    %[b,bint]=regress(yor(1:(29+iy)),Xnat(1:(29+iy),:));
    %Bnat(iy,1)=b(2);
    %Bnat(iy,2:3)=bint(2,:);
    bootdat = [yor(1:(29+iy)),Xnat(1:(29+iy),:)];
    aux=bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
    Bnat_boot(iy,:,:) = aux;
end

indices{1}=[1:size(Ball_boot,2)];
boundary=cell(1,1);boundary{1}={'prc5','prc95'}; 
col_all_spread=[237 220 202]/255;
col_all_line=[130 90 38]/255;
col_nat_spread=[181 209 206]/255;
col_nat_line=[65 116 110]/255;

figure; hold on
plot(2000:2021,prctile(Ball_boot(:,:,2),[50],2),'-','LineWidth',2,'Color',col_all_line)
plot(2000:2021,prctile(Bnat_boot(:,:,2),[50],2),'-','LineWidth',2,'Color',col_nat_line)
legend('CMIP6-ALL','CMIP6-NAT','Location','Best','AutoUpdate','off')
legend boxoff
drawSpread(Bnat_boot(:,:,2),'xvalues',(2000:2021)','indexesg',indices,'colorsg',col_nat_spread,'lines','no','boundary',boundary,'alphasg',[0.5 0.5])
drawSpread(Ball_boot(:,:,2),'xvalues',(2000:2021)','indexesg',indices,'colorsg',col_all_spread,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]);
xlabel('Last year of the serie','FontSize',18) 
ylabel('Scaling factor','FontSize',18) 
set(gca,'FontSize',18,'xlim',[2000,2021])
gridxy([],[0],'Color',[0 0 0]/255,'Linestyle','-','LineWidth',2);
ylim([-2.5 2.5])
file=[dir_out,'scaling_factor_1member.eps']
% print( gcf, '-depsc2', file ,'-painters')
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')

%% boxplots impacts

pred_all_12_reshaped=reshape(pred_all_12,[size(pred_all_12,1),size(pred_all_12,2)*size(pred_all_12,3)]);
pred_nat_12_reshaped=reshape(pred_nat_12,[size(pred_nat_12,1),size(pred_nat_12,2)*size(pred_nat_12,3)]);
size(pred_all_12_reshaped)

for ib=1:size(pred_all_12_reshaped,2)
    %prctile(nanmean(exp(pred_all_12(:,:)),1),[2.5 50 97.5],2)
    %prctile(nanmean(exp(pred_nat_12(:,:)),1),[2.5 50 97.5],2)
    impacts_50(ib)=100*(nanmean(exp(pred_all_12_reshaped(:,ib)),1)-nanmean(exp(pred_nat_12_reshaped(:,ib)),1))./nanmean(exp(pred_nat_12_reshaped(:,ib)),1);
    impacts_fh(ib)=100*(nanmean(exp(pred_all_12_reshaped(1:25,ib)),1)-nanmean(exp(pred_nat_12_reshaped(1:25,ib)),1))./nanmean(exp(pred_nat_12_reshaped(1:25,ib)),1);
    impacts_sh(ib)=100*(nanmean(exp(pred_all_12_reshaped(26:51,ib)),1)-nanmean(exp(pred_nat_12_reshaped(26:51,ib)),1))./nanmean(exp(pred_nat_12_reshaped(26:51,ib)),1);
end


round(prctile(impacts_fh,[2.5 50 97.5],2))
round(prctile(impacts_sh,[2.5 50 97.5],2))
round(prctile(impacts_50,[2.5 50 97.5],2))




%%  boxplot 3 periods
dataplot=zeros(1,size(impacts_50,2),3)*NaN;

dataplot(1,1:size(impacts_fh,2),1)=impacts_fh';
dataplot(1,1:size(impacts_sh,2),2)=impacts_sh';
dataplot(1,1:size(impacts_50,2),3)=impacts_50';

xlabel_periods(1,:)=['1971 - 1995'];
xlabel_periods(2,:)=['1996 - 2021'];
xlabel_periods(3,:)=['1971 - 2021'];

%c=gray(3);
%c=c(2,:);
c=[0 102 204]/255;
hFig=figure;aboxplot3(dataplot,'colormap',c)
hold on
%gridxy([size(impact_20y,1)+1],[],'Color',[100 100 100]/255,'Linestyle','-');
gridxy([],[-100:100:1200],'Color',[121 121 121]/255,'Linestyle',':');
gridxy([],[0],'Color',[0 0 0],'Linestyle','-','LineWidth',2);
ylabel('Impact [%]','FontSize',font_size);
xlabel('periods','FontSize',font_size);
set(gca,'XTick',1:size(xlabel_periods,1));
set(gca,'XTickLabel',xlabel_periods,'Fontsize',font_size);
htick=xticklabel_rotate([],45,[]);
set(gca,'FontSize',font_size)

nomeout=[dir_out,'figure_3c_boxplot_3.eps'];
print( gcf, '-depsc2', nomeout ,'-painters')


%% boxplots impacts 1 member

pred_all_1_reshaped=reshape(pred_all_1,[size(pred_all_1,1),size(pred_all_1,2)*size(pred_all_1,3)]);
pred_nat_1_reshaped=reshape(pred_nat_1,[size(pred_nat_1,1),size(pred_nat_1,2)*size(pred_nat_1,3)]);


for ib=1:size(pred_all_1_reshaped,2)
    impacts_50(ib)=100*(nanmean(exp(pred_all_1_reshaped(:,ib)),1)-nanmean(exp(pred_nat_1_reshaped(:,ib)),1))./nanmean(exp(pred_nat_1_reshaped(:,ib)),1);
    impacts_fh(ib)=100*(nanmean(exp(pred_all_1_reshaped(1:25,ib)),1)-nanmean(exp(pred_nat_1_reshaped(1:25,ib)),1))./nanmean(exp(pred_nat_1_reshaped(1:25,ib)),1);
    impacts_sh(ib)=100*(nanmean(exp(pred_all_1_reshaped(26:51,ib)),1)-nanmean(exp(pred_nat_1_reshaped(26:51,ib)),1))./nanmean(exp(pred_nat_1_reshaped(26:51,ib)),1);
end


round(prctile(impacts_fh,[2.5 50 97.5],2))
round(prctile(impacts_sh,[2.5 50 97.5],2))
round(prctile(impacts_50,[2.5 50 97.5],2))




%%  boxplot 3 periods 1 member
dataplot=zeros(1,size(impacts_50,2),3)*NaN;

dataplot(1,1:size(impacts_fh,2),1)=impacts_fh';
dataplot(1,1:size(impacts_sh,2),2)=impacts_sh';
dataplot(1,1:size(impacts_50,2),3)=impacts_50';

xlabel_periods(1,:)=['1971 - 1995'];
xlabel_periods(2,:)=['1996 - 2021'];
xlabel_periods(3,:)=['1971 - 2021'];

%c=gray(3);
%c=c(2,:);
c=[0 102 204]/255;
hFig=figure;aboxplot3(dataplot,'colormap',c)
hold on
%gridxy([size(impact_20y,1)+1],[],'Color',[100 100 100]/255,'Linestyle','-');
gridxy([],[-100:100:1200],'Color',[121 121 121]/255,'Linestyle',':');
gridxy([],[0],'Color',[0 0 0],'Linestyle','-','LineWidth',2);
ylabel('Impact [%]','FontSize',font_size);
xlabel('periods','FontSize',font_size);
set(gca,'XTick',1:size(xlabel_periods,1));
set(gca,'XTickLabel',xlabel_periods,'Fontsize',font_size);
htick=xticklabel_rotate([],45,[]);
set(gca,'FontSize',font_size)

nomeout=[dir_out,'figure_3c_boxplot_3_1member.eps'];
print( gcf, '-depsc2', nomeout ,'-painters')

%% PRECIPITATION

clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years_study=1971:2021;
years_nat_all=1850:2021;
years_nat=1950:2021;
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
font_size=12;
base_period=1960:1982;
best_start_tx=4;
best_stop_tx=10;
alpha = 0.05;
NB=10000;
n=1; %linear trend

%% load obs climate 


%climatic variables
filename = [dir_data,'nclimgrid/nclimgrid.mat'];
load(filename);

PREC_obs = zeros(length(years_climate),1)*NaN;
PREC_obs_year = zeros(length(years_climate),1)*NaN;
for iyear=1:length(years_climate) 
  i1 = (iyear - 1) * 12 + best_start_tx;
  i2 = (iyear - 1) * 12 + best_stop_tx;
  PREC_obs(iyear) = sum(nclimgrid(2,i1:i2));
  i1 = (iyear - 1) * 12 + 1;
  i2 = (iyear - 1) * 12 + 12;
  PREC_obs_year(iyear) = sum(nclimgrid(2,i1:i2));
end
%figure;plot(TSMAX_obs)
PREC_obs_anom=scale_base_period(PREC_obs,base_period,years_climate);
PREC_obs_year_anom=scale_base_period(PREC_obs_year,base_period,years_climate);
[C,IA,IB] = intersect(years_study,years_climate);
PREC_obs_year_anom=PREC_obs_year_anom(IB);
PREC_obs_anom=PREC_obs_anom(IB);
PREC_obs=PREC_obs(IB);

mean(PREC_obs(26:51))-mean(PREC_obs(1:25))
mean(PREC_obs_anom(26:51))-mean(PREC_obs_anom(1:25))
mean(PREC_obs_year(26:51))-mean(PREC_obs_year(1:25))

mean(PREC_obs_year)
mean(PREC_obs)

mean(PREC_obs)/mean(PREC_obs_year)
%% load data hist-nat vs all
nomefile=[dir_data 'gcms/GCMs_attribution_prec.mat'];
load(nomefile)

%iok2=sum(TSMAX_ALL,1)
%length(find(isnan(iok2)))
aux=sum(PREC_NAT,1);
iok=(find(~isnan(aux)));
aux=cell(1,length(iok));
for ii = 1:length(iok)
  aux{ii} = model_com{iok(ii)};
end
model_com=aux;
model_member_num=model_member_num(iok,:);
PREC_ALL=PREC_ALL(:,iok);
PREC_NAT=PREC_NAT(:,iok);
PREC_ALL_ANOM=PREC_ALL_ANOM(:,iok);
PREC_NAT_ANOM=PREC_NAT_ANOM(:,iok);

for im=1:model_member_num(end,1)
    num_member=find(model_member_num(:,1)==im);
    PREC_ALL_12(:,im)=nansum(PREC_ALL(:,num_member),2)/length(num_member);
    PREC_NAT_12(:,im)=nansum(PREC_NAT(:,num_member),2)/length(num_member);
    PREC_ALL_12_ANOM(:,im)=nansum(PREC_ALL_ANOM(:,num_member),2)/length(num_member);
    PREC_NAT_12_ANOM(:,im)=nansum(PREC_NAT_ANOM(:,num_member),2)/length(num_member);
end


PREC_NAT_12_ANOM=PREC_NAT_12_ANOM(IB,:);
PREC_ALL_12_ANOM=PREC_ALL_12_ANOM(IB,:);
PREC_NAT_12=PREC_NAT_12(IB,:);
PREC_ALL_12=PREC_ALL_12(IB,:);

figure; plot(movmean(PREC_NAT_12_ANOM,21))
        

col_all_spread=[237 220 202]/255;
col_all_line=[130 90 38]/255;
col_nat_spread=[181 209 206]/255;
col_nat_line=[65 116 110]/255;

figure; hold on
plot(years_study,PREC_obs(:),'k:o','LineWidth',1)
plot(years_study,PREC_NAT_12,'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,PREC_ALL_12,'-','LineWidth',0.5,'Color',col_all_spread)



%% PREC AVERAGED MEMBERS FOR EACH MODEL


trend_obs=movmean(PREC_obs_anom,21);
trend_nat=movmean(nanmean(PREC_NAT_12_ANOM,2),21); 
trend_all=movmean(nanmean(PREC_ALL_12_ANOM,2),21); 


indices{1}=[1:size(PREC_ALL_12,2)];
boundary=cell(1,1);boundary{1}={'prc2.5','prc97.5'}; 

figure; hold on
plot(years_study,(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
legend('CMIP6-ALL','CMIP6-NAT','Observed','Location','Best','AutoUpdate','off')
legend boxoff
plot(years_study,(PREC_NAT_12_ANOM),'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,(PREC_ALL_12_ANOM),'-','LineWidth',0.5,'Color',col_all_spread)
plot(years_study,(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(years_study,PREC_obs_anom,'k:o','LineWidth',1)
xlabel('years','FontSize',font_size) 
ylabel('Precipitation [mm]','FontSize',font_size) 
set(gca,'FontSize',font_size,'xlim',[years_study(1)-1 years_study(end)+1])
%gridxy([],[-2.5:0.5:2.5],'Color',[121 121 121]/255,'Linestyle','-');
%ylim([-2.5 2.5])
file=[dir_out,'figure_2a_prec_all_members.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')





(mean(PREC_obs_anom(26:51)))-(mean(PREC_obs_anom(1:25)))
prctile((mean(PREC_NAT_12_ANOM(26:51,:),1)-mean(PREC_NAT_12_ANOM(1:25,:),1)),[2.5 50 97.5],2)
prctile((mean(PREC_ALL_12_ANOM(26:51,:),1)-mean(PREC_ALL_12_ANOM(1:25,:),1)),[2.5 50 97.5],2)

%% VPD
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts_def/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years_study=1971:2021;
years_nat_all=1850:2021;
years_nat=1950:2021;
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
font_size=12;
base_period=1960:1982;
best_start_tx=4;
best_stop_tx=10;
alpha = 0.05;
NB=10000;
n=1; %linear trend

%% load obs
%climatic variables prism
filename = [dir_data,'prism/vpd_prism.mat']; %period 1960:2021
load(filename);
%vpd_prism=vpd_prism((11*12+1):end); %period 1971:2021


VPD_obs = zeros(length(years_climate),1)*NaN;
for iyear=1:length(years_climate)
    i1 = (iyear - 1) * 12 + best_start_tx;
    i2 = (iyear - 1) * 12 + best_stop_tx;
    VPD_obs(iyear) = mean(vpd_prism(i1:i2));
end
%figure;plot(TSMAX_obs)
VPD_obs_anom=scale_base_period(VPD_obs,base_period,years_climate);
%figure;plot(VPD_obs_anom)
figure;plot(VPD_obs*10)
[C,IA,IB] = intersect(years_study,years_climate);
VPD_obs_anom=VPD_obs_anom(IB)*10; %hPa
VPD_obs=VPD_obs(IB);
mean(VPD_obs(26:51))-mean(VPD_obs(1:25))



%% load data hist-nat vs all
nomefile=[dir_data 'gcms/GCMs_attribution_vpd.mat'];
load(nomefile)
[C,IA,IB]=intersect(model_nat,model_all);
model_nat=model_nat(1,IA);
VPD_NAT=VPD_NAT(:,IA);
model_all=model_all(1,IB);
VPD_ALL=VPD_ALL(:,IB);



aux=sum(VPD_NAT,1);
iok=(find(~isnan(aux)));
aux=cell(1,length(iok));
for ii = 1:length(iok)
  aux{ii} = model_nat{iok(ii)};
end
model_nat=aux;
model_all=aux;
VPD_ALL=VPD_ALL(:,iok);
VPD_NAT=VPD_NAT(:,iok);
VPD_ALL_ANOM=VPD_ALL_ANOM(:,iok);
VPD_NAT_ANOM=VPD_NAT_ANOM(:,iok);

aux2=split(model_nat,";");
single_model_ok=unique(aux2(1,:,1));


figure;plot(VPD_ALL)
figure;plot(VPD_NAT)
for im=1:length(single_model_ok)
    single_model_ok(im)
    num_member=find(strcmp(aux2(1,:,1),single_model_ok(im)));
    VPD_ALL_12(:,im)=sum(VPD_ALL(:,num_member),2)/length(num_member);
    VPD_NAT_12(:,im)=sum(VPD_NAT(:,num_member),2)/length(num_member);
    VPD_ALL_12_ANOM(:,im)=sum(VPD_ALL_ANOM(:,num_member),2)/length(num_member);
    VPD_NAT_12_ANOM(:,im)=sum(VPD_NAT_ANOM(:,num_member),2)/length(num_member);
end

[C,IA,IB] = intersect(years_study,years_climate);
VPD_NAT_12_ANOM=VPD_NAT_12_ANOM(IB,:)*10; %hPa
VPD_ALL_12_ANOM=VPD_ALL_12_ANOM(IB,:)*10; %hPa
VPD_NAT_12=VPD_NAT_12(IB,:);
VPD_ALL_12=VPD_ALL_12(IB,:);

col_all_spread=[237 220 202]/255;
col_all_line=[130 90 38]/255;
col_nat_spread=[181 209 206]/255;
col_nat_line=[65 116 110]/255;

figure; hold on
plot(years_study,VPD_obs_anom(:),'k:o','LineWidth',1)
plot(years_study,VPD_NAT_12_ANOM,'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,VPD_ALL_12_ANOM,'-','LineWidth',0.5,'Color',col_all_spread)

figure; hold on
plot(years_study,VPD_obs(:),'k:o','LineWidth',1)
plot(years_study,VPD_NAT_12,'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,VPD_ALL_12,'-','LineWidth',0.5,'Color',col_all_spread)

%% VPD AVERAGED MEMBERS FOR EACH MODEL


trend_obs=movmean(VPD_obs_anom,21);
trend_nat=movmean(nanmean(VPD_NAT_12_ANOM,2),21); 
trend_all=movmean(nanmean(VPD_ALL_12_ANOM,2),21); 

indices{1}=[1:size(VPD_ALL_12,2)];
boundary=cell(1,1);boundary{1}={'prc2.5','prc97.5'}; 

figure; hold on
plot(years_study,(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
legend('CMIP6-ALL','CMIP6-NAT','Observed','Location','Best','AutoUpdate','off')
legend boxoff
plot(years_study,(VPD_NAT_12_ANOM),'-','LineWidth',0.5,'Color',col_nat_spread)
plot(years_study,(VPD_ALL_12_ANOM),'-','LineWidth',0.5,'Color',col_all_spread)
plot(years_study,(trend_all),'-','LineWidth',2,'Color',col_all_line)
plot(years_study,(trend_nat),'-','LineWidth',2,'Color',col_nat_line)
plot(years_study,(trend_obs),'k-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(years_study,VPD_obs_anom,'k:o','LineWidth',1)
xlabel('years','FontSize',font_size) 
ylabel('VPD [hPa]','FontSize',font_size) 
set(gca,'FontSize',font_size,'xlim',[years_study(1)-1 years_study(end)+1])
%gridxy([],[-2.5:0.5:2.5],'Color',[121 121 121]/255,'Linestyle','-');
%ylim([-2.5 2.5])
file=[dir_out,'figure_2a_vpd_all_members.eps']
set(gcf,'PaperType','A4')
print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')


((mean(VPD_obs(26:51)))-(mean(VPD_obs(1:25))))*10
((mean(VPD_obs_anom(26:51)))-(mean(VPD_obs_anom(1:25))))
prctile((mean(VPD_NAT_12_ANOM(26:51,:),1)-mean(VPD_NAT_12_ANOM(1:25,:),1)),[2.5 50 97.5],2)
prctile((mean(VPD_ALL_12_ANOM(26:51,:),1)-mean(VPD_ALL_12_ANOM(1:25,:),1)),[2.5 50 97.5],2)
