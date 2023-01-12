%% plot figure 1b, 1c, test best model
clear all,clc,close all
addpath '~/Dropbox/estcena/scripts/fires_california/scripts/misc/'
savepath
workPath=[pwd];cd([workPath])

%% misc
years_climate=1960:2021;
years=1971:2021;
base_period=1960:1982;
dir_data='~/Dropbox/estcena/scripts/fires_california/data_def/';
dir_out='~/Dropbox/estcena/scripts/fires_california/paper/figs_temp/';
alpha = 0.05;
NB=10000;
best_start_tx=4;
best_stop_tx=10;
font_size=12;

%% load data

%frap 
namefile = [dir_data,'fires/frap_forest_sierra_ncoast_year.mat'];
load(namefile,'FIRE')

%climatic variables nclimgrid
filename = [dir_data,'nclimgrid/nclimgrid.mat']; %period 1960:2021
load(filename);
nclimgrid=nclimgrid(:,(11*12+1):end); %period 1971:2021

%climatic variables prism
filename = [dir_data,'prism/vpd_prism.mat']; %period 1960:2021
load(filename);
vpd_prism=vpd_prism((11*12+1):end); %period 1971:2021
filename = [dir_data,'prism/tmax_prism.mat']; %period 1960:2021
load(filename);
tmax_prism=tmax_prism((11*12+1):end); %period 1971:2021

figure;plot(tmax_prism,nclimgrid(1,:),'o')
corrcoef(tmax_prism,nclimgrid(1,:))

VPD_3_10 = zeros(length(years),1)*NaN;
VPD_4_10 = zeros(length(years),1)*NaN;
VPD_5_10 = zeros(length(years),1)*NaN;
VPD_5_9 = zeros(length(years),1)*NaN;

TSMAX_3_10 = zeros(length(years),1)*NaN;
TSMAX_4_10 = zeros(length(years),1)*NaN;
TSMAX_5_10 = zeros(length(years),1)*NaN;
TSMAX_5_9 = zeros(length(years),1)*NaN;

PREC_3_10 = zeros(length(years),1)*NaN;
PREC_4_10 = zeros(length(years),1)*NaN;
PREC_5_10 = zeros(length(years),1)*NaN;
PREC_5_9 = zeros(length(years),1)*NaN;
PREC_1_12 = zeros(length(years),1)*NaN;

for iyear=1:length(years) 
  i1 = (iyear - 1) * 12 + best_start_tx;
  i2 = (iyear - 1) * 12 + best_stop_tx;
  TSMAX_4_10(iyear) = mean(nclimgrid(1,i1:i2));
  PREC_4_10(iyear) = mean(nclimgrid(2,i1:i2));
  VPD_4_10(iyear) = mean(vpd_prism(i1:i2));
  i1 = (iyear - 1) * 12 + 3;
  i2 = (iyear - 1) * 12 + 10;
  TSMAX_3_10(iyear) = mean(nclimgrid(1,i1:i2));
  PREC_3_10(iyear) = mean(nclimgrid(2,i1:i2));
  VPD_3_10(iyear) = mean(vpd_prism(i1:i2));
  i3 = (iyear - 1) * 12 + 5;
  i4 = (iyear - 1) * 12 + 10;
  TSMAX_5_10(iyear) = mean(nclimgrid(1,i3:i4));  
  PREC_5_10(iyear) = mean(nclimgrid(2,i3:i4)); 
  VPD_5_10(iyear) = mean(vpd_prism(i3:i4));
  i3 = (iyear - 1) * 12 + 5;
  i4 = (iyear - 1) * 12 + 9;
  TSMAX_5_9(iyear) = mean(nclimgrid(1,i3:i4)); 
  PREC_5_9(iyear) = mean(nclimgrid(2,i3:i4)); 
  VPD_5_9(iyear) = mean(vpd_prism(i3:i4));
  i3 = (iyear - 1) * 12 + 1;
  i4 = (iyear - 1) * 12 + 12;
  PREC_1_12(iyear) = sum(nclimgrid(2,i3:i4));  
end

mean(PREC_1_12(26:51))-mean(PREC_1_12(1:25))
TSMAX_4_10_orig=TSMAX_4_10;
PREC_4_10_orig=PREC_4_10;

%figure;plot(TSMAX)
TSMAX_3_10=scale_base_period(TSMAX_3_10,base_period,years);
TSMAX_4_10=scale_base_period(TSMAX_4_10,base_period,years);
TSMAX_5_10=scale_base_period(TSMAX_5_10,base_period,years);
TSMAX_5_9=scale_base_period(TSMAX_5_9,base_period,years);
%figure;plot(TSMAX)
PREC_3_10=scale_base_period(PREC_3_10,base_period,years);
PREC_4_10=scale_base_period(PREC_4_10,base_period,years);
PREC_5_10=scale_base_period(PREC_5_10,base_period,years);
PREC_5_9=scale_base_period(PREC_5_9,base_period,years);

%% correlations
clc
[rho,sig]=corr_boot(log(FIRE),TSMAX_5_9)
[rho,sig]=corr_boot(log(FIRE),TSMAX_5_10)
[rho,sig]=corr_boot(log(FIRE),TSMAX_4_10)
[rho,sig]=corr_boot(log(FIRE),TSMAX_3_10)
clc
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(TSMAX_5_9))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(TSMAX_5_10))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(TSMAX_4_10))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(TSMAX_3_10))
clc
[rho,sig]=corr_boot(log(FIRE),VPD_5_9)
[rho,sig]=corr_boot(log(FIRE),VPD_5_10)
[rho,sig]=corr_boot(log(FIRE),VPD_4_10)
[rho,sig]=corr_boot(log(FIRE),VPD_3_10)
clc
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(VPD_5_9))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(VPD_5_10))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(VPD_4_10))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(VPD_3_10))
clc
[rho,sig]=corr_boot((log(FIRE)),(PREC_5_9))
[rho,sig]=corr_boot((log(FIRE)),(PREC_5_10))
[rho,sig]=corr_boot((log(FIRE)),(PREC_4_10))
[rho,sig]=corr_boot((log(FIRE)),(PREC_3_10))
clc
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(PREC_5_9))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(PREC_5_10))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(PREC_4_10))
[rho,sig]=corr_boot(detrend(log(FIRE)),detrend(PREC_3_10))
clc
[rho,sig]=corr_boot(PREC_5_9,VPD_5_9)
[rho,sig]=corr_boot(PREC_5_10,VPD_5_10)
[rho,sig]=corr_boot(PREC_4_10,VPD_4_10)
[rho,sig]=corr_boot(PREC_3_10,VPD_3_10)
clc
[rho,sig]=corr_boot(detrend(PREC_5_9),detrend(VPD_5_9))
[rho,sig]=corr_boot(detrend(PREC_5_10),detrend(VPD_5_10))
[rho,sig]=corr_boot(detrend(PREC_4_10),detrend(VPD_4_10))
[rho,sig]=corr_boot(detrend(PREC_3_10),detrend(VPD_3_10))
clc
[rho,sig]=corr_boot(TSMAX_5_9,VPD_5_9)
[rho,sig]=corr_boot(TSMAX_5_10,VPD_5_10)
[rho,sig]=corr_boot(TSMAX_4_10,VPD_4_10)
[rho,sig]=corr_boot(TSMAX_3_10,VPD_3_10)
clc
[rho,sig]=corr_boot(detrend(TSMAX_5_9),detrend(VPD_5_9))
[rho,sig]=corr_boot(detrend(TSMAX_5_10),detrend(VPD_5_10))
[rho,sig]=corr_boot(detrend(TSMAX_4_10),detrend(VPD_4_10))
[rho,sig]=corr_boot(detrend(TSMAX_3_10),detrend(VPD_3_10))
clc
[rho,sig]=corr_boot(TSMAX_5_9,PREC_5_9)
[rho,sig]=corr_boot(TSMAX_5_10,PREC_5_10)
[rho,sig]=corr_boot(TSMAX_4_10,PREC_4_10)
[rho,sig]=corr_boot(TSMAX_4_10_orig,PREC_4_10_orig)
[rho,sig]=corr_boot(TSMAX_3_10,PREC_3_10)
clc
[rho,sig]=corr_boot(detrend(TSMAX_5_9),detrend(PREC_5_9))
[rho,sig]=corr_boot(detrend(TSMAX_5_10),detrend(PREC_5_10))
[rho,sig]=corr_boot(detrend(TSMAX_4_10),detrend(PREC_4_10))
[rho,sig]=corr_boot(detrend(TSMAX_3_10),detrend(PREC_3_10))
clc
[rho,sig] = partialcorr(TSMAX_4_10,log(FIRE),PREC_4_10)


%% regression model tsmax prec 1971-2021
yor=log(FIRE);
Xorig1 = [ones(size(TSMAX_4_10,1),1) TSMAX_4_10 PREC_4_10]; %osservati
[B1,BINT,res0,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(res0),nanvar((yor)))
bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootCI1=prctile(bootb1,[2.5 97.5]);
format longG
[B1(2) bootCI1(1,2)'  bootCI1(2,2)']
[B1(3) bootCI1(1,3)'  bootCI1(2,3)']


%% plot Fig. 1a
x1=years;
x2=x1;
y1=FIRE;
y2=TSMAX_4_10_orig;


%figure;
figure1 = figure;
colormap(white);

% Create axes
axes1 = axes('Parent',figure1);
grayColor = [.3 .3 .3];
%hl1 = bar(x1,y1,'w');
hl1 = line(x1,y1,'LineStyle','-','Marker','o','MarkerSize',10, 'MarkerFaceColor','r','Color','r');
%line(x1(1:25),repmat(mean(y1(1:25)),1,25),'LineStyle',':','Color',grayColor)
set(gca,'XTick',[])
xlim([x1(1)-1 x1(end)+1]);
%ylim([0 11]);
%set(hl1, 'YTick', [1:1:10])
colormap white
legend('BA','Location','NorthWest')
ax1 = gca;
set(ax1, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     )
set(ax1, 'YScale', 'log')
%set(gca,'FontSize',font_size)
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
hl2 = line(x2,y2,'LineStyle','-','Marker','o','MarkerSize',10, 'MarkerFaceColor','w','Color','k','Parent',ax2);
%line(x1(26:60),repmat(mean(y1(26:50)),1,25),'LineStyle',':','Color',grayColor,'Parent',ax2)
%ylim([21 27]);
xlabel('Years','FontSize',font_size)
set(get(ax1,'Ylabel'),'String','Burned Area (km^2)','color','k','FontSize',font_size)
set(get(ax2,'Ylabel'),'String','TS_{max} (??C)','Color','k','FontSize',font_size)
%h=legend('TS_{max}','Location','NorthEast');
h=legend('TS_{max}','Location','SouthEast');
set(h, 'Color', 'none')
%set(gca,'FontSize',18)
set(ax2,'TickDir','Out')
ylimits1 = get(ax1,'YLim');
ylimits2 = get(ax2,'YLim');
xlim([x1(1)-1 x1(end)+1]);
%set(ax1,'YTick',[0:1000:1000])
set(ax1,'FontSize',font_size,'TickDir','out','XTick',[],'YMinorTick','on',...
    'YScale','log','YTick',[10 25 50 100  200 500  1000  2500  5000 9500]);
set(ax2,'YTick',[ylimits2(1):1:ylimits2(2)])
set(gca,'FontSize',font_size)
set(gcf,'PaperType','A4')
file=[dir_out,'figure_1a_temp.pdf']
print( gcf, '-dpdf', file ,'-painters')
%file=[dir_out,'fig1a.eps']
%print( gcf, '-depsc2', file ,'-painters')

mean(TSMAX_4_10_orig(1:25))
mean(TSMAX_4_10_orig(26:51))

mean(FIRE(1:25))
mean(FIRE(26:51))

mean(FIRE(26:51))/mean(FIRE(1:25))

mean(FIRE(end-19:end))



%% regression model tsmax 1971-2021
yor=log(FIRE);
% regression only TSMAX
Xorig1 = [ones(size(TSMAX_4_10,1),1) TSMAX_4_10]; %osservati
[B1,BINT,res0,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootCI1=prctile(bootb1,[2.5 97.5]);
format longG
[B1(2) bootCI1(1,2)'  bootCI1(2,2)']
100*(exp([B1(2) bootCI1(1,2)'  bootCI1(2,2)'])-1)

Y = (Xorig1*B1);

% testing residuals tsmax model

%The Durbin-Watson statistic [3] is the autocorrelation measure most frequently reported in econometric analyses. One reason is that it is easy to compute. For the M0 model:
%[3] Durbin, J. and G.S. Watson. "Testing for Serial Correlation in Least Squares Regression." Biometrika. Vol. 37, 1950, pp. 409?428.
p1 = dwtest(res0,Xorig1)
%0.8943 The p-value for a null of no first-order autocorrelation is well above the standard 5% critical value.
%[4] Engle, Robert. F. ?Autoregressive Conditional Heteroscedasticity with Estimates of the Variance of United Kingdom Inflation.? Econometrica 50 (July 1982): 987?1007. https://doi.org/10.2307/1912773.
[hARCH0,pARCH0] = archtest(res0)
%One-sample Kolmogorov-Smirnov test
[H,P] = kstest(res0)

%% out-of-sample figure 1b

oos2=[]; 
p=10; %leave-ten-years-out
numanni=length(years);
Xorig = [ones(size(TSMAX_4_10,1),1)  TSMAX_4_10 ]; %osservati
for i=1:p:(numanni-1);
    if i==41
        itest=[i:i+(p)];
    else
        itest=[i:i+(p-1)];
    end
    itrain1=setdiff([1:numanni],itest);
    X=Xorig(itrain1,:); 
    y=yor(itrain1);
    [B,BINT,R,RINT,STATS] = regress(y,X,alpha);
    sigma=std(R);
    clear X
    X=Xorig(itest,:); 
    oos1=[]; 
    for iboot=1:NB
        noise=normrnd(0,sigma,length(itest),1);
        Y = (X*B)+noise;
        oos1=[oos1,Y];
    end   
    oos2=[oos2;oos1];
end    

for i=1:numanni
    %b2(i,:)=[min(exp(oos2(i,:))-nf(i)),max(exp(oos2(i,:))-nf(i))];
    Ymed(i)=prctile(exp(oos2(i,:)),50);
    b1(i,:)=[prctile(exp(oos2(i,:)),2.5),prctile(exp(oos2(i,:)),97.5)];
    
end 
figure1 = figure('PaperSize',[20.98 29.68]);
%ylim([-4 4])
%xlim([-4 4])
% Create axes
axes1 = axes('Parent',figure1);
hold('all');
neg=b1(:,1);
pos=b1(:,2);
for ip=1:length(Ymed)
    ly=neg(ip);
    uy=pos(ip);
    %l1=line([log(BAS(ip)) log(BAS(ip))],[ly uy],'color',[128 128 128]/255,'LineWidth',0.5);
    l1=line([exp(yor(ip)) exp(yor(ip))],[ly uy],'color',[128 128 128]/255,'LineWidth',0.5);
    %l2=line([TMEAN_ECO_YM(ip)-0.1 TMEAN_ECO_YM(ip)+0.1],[ly ly],'color',[128 128 128]/255,'LineWidth',0.5);
    %l3=line([TMEAN_ECO_YM(ip)-0.1 TMEAN_ECO_YM(ip)+0.1],[uy uy],'color',[128 128 128]/255,'LineWidth',0.5);
end
ind_years=years*NaN;
ind_years(find(years>=1971 & years<=1980))=1;
ind_years(find(years>=1980 & years<=1990))=2;
ind_years(find(years>=1991 & years<=2000))=3;
ind_years(find(years>=2001 & years<=2010))=4;
ind_years(find(years>=2011 & years<=2021))=5;
scatter1 =scatter(exp(yor),(Ymed),150,ind_years,'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%caxis([1971 2020])
maxdensity=5;
cmap = flipud(autumn(maxdensity));
colormap(cmap)

h=colorbar

%title(h,'Year','FontSize',font_size)
%YTicknome={ '1971-1980', '1981-1990','1991-2000','2001-2010','2011-2020'};
YTicknome={ '1971', '1980','1990','2000','2010','2020'};
% set(h, 'YTick', [1:maxdensity])
% set(h,'YTickLabel',{[ 'NaN '], [ '<1 '] , ...
%     [ '1-15'] , ...
%     [ '15-30'] , [ '30-45'] ,...
%     [ '>45'] }) 
% cbarHandle = colorbar('YTick',...
% [1+0.5*(maxdensity-1)/maxdensity:(maxdensity-1)/maxdensity:maxdensity],...
% 'YTickLabel',YTicknome, 'YLim', [1 maxdensity]);
cbarHandle = colorbar('YTick',...
[1:(maxdensity-1)/maxdensity:maxdensity],...
'YTickLabel',YTicknome, 'YLim', [1 maxdensity]);

scatter1.MarkerFaceAlpha = 0.75;

line([exp(1.8) exp(10.3)],[exp(1.8) exp(10.3)],'color','k')
xlim([exp(1.8) exp(10.3)])
ylim([exp(1.8) exp(10.3)])
xlabel('Observed BA [km^2]','FontSize',font_size)
ylabel('Predicted BA [km^2]','FontSize',font_size)
set(gca,'FontSize',font_size)
set(gcf,'PaperType','A4')
file=[dir_out,'figure_1b_temp.pdf']
print( gcf, '-dpdf', file ,'-painters')
%file=[dir_out,'figure1d_temp.eps']
%print( gcf, '-depsc2', file ,'-painters')

[r p]=corrcoef((yor),log(Ymed),'rows','complete')
[r p]=corrcoef((yor),TSMAX_4_10,'rows','complete')

%% stability

yor=log(FIRE);
Xorig1 = [ones(size(TSMAX_4_10,1),1) (TSMAX_4_10)]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
bootdat1 = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
bootCI1=prctile(bootb1,[2.5 97.5])
B1

%detrend
yor=log(FIRE);
Xorig1 = [ones(size(TSMAX_4_10,1),1) detrend(TSMAX_4_10)]; %osservati
[B_detrend,BINT,R,RINT,STATS1] = regress(detrend(yor),Xorig1,alpha);
STATS1
bootdat1 = [detrend(yor),Xorig1];
bootb_detrend = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
bootCI_detrend=prctile(bootb_detrend,[2.5 97.5])
B_detrend

%first-half period
Xorig1 = [ones(size(TSMAX_4_10(1:25),1),1) (TSMAX_4_10(1:25))]; %osservati
[B11] = regress((yor(1:25)),Xorig1,alpha);
STATS1
bootdat1 = [(yor(1:25)),Xorig1];
bootb11 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
bootCI11=prctile(bootb11,[2.5 97.5])

%second-half period
Xorig1 = [ones(size(TSMAX_4_10(26:51),1),1) (TSMAX_4_10(26:51))]; %osservati
[B21] = regress((yor(26:51)),Xorig1,alpha);
bootdat1 = [(yor(26:51)),Xorig1];
bootb21 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
bootCI21=prctile(bootb21,[2.5 97.5])

%25 coldest years
nextremes=25;
[~,i_sorted]=sort(TSMAX_4_10);
Xorig1 = [ones(size(TSMAX_4_10(1:nextremes),1),1) (TSMAX_4_10(i_sorted(1:nextremes)))]; %osservati
[B_smallest] = regress((yor(i_sorted(1:nextremes))),Xorig1,alpha);
STATS1
bootdat1 = [(yor(i_sorted(1:nextremes))),Xorig1];
bootb_coldest = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
bootCI_smallest=prctile(bootb_coldest,[2.5 97.5])

%25 hottest years
[~,i_sorted]=sort(TSMAX_4_10,'descend');
Xorig1 = [ones(size(TSMAX_4_10(1:nextremes),1),1) (TSMAX_4_10(i_sorted(1:nextremes)))]; %osservati
[B_largest] = regress((yor(i_sorted(1:nextremes))),Xorig1,alpha);
STATS1
bootdat1 = [(yor(i_sorted(1:nextremes))),Xorig1];
bootb11 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);
bootCI_largest=prctile(bootb11,[2.5 97.5])

xlabel_periods(1,:)=['1971 - 2021'];
xlabel_periods(2,:)=['detrended  '];
xlabel_periods(3,:)=['1971 - 1995'];
xlabel_periods(4,:)=['1996 - 2021'];
xlabel_periods(5,:)=['25 coldest '];
xlabel_periods(6,:)=['25 hottest '];

for ip=1:5
    
    ip1=((ip-1)*5)+1;
    if ip==5
        ip2=ip1+30;
    else
        ip2=ip1+29;
    end
    xlabel_periods(ip+6,:)=[num2str(years(ip1)),' - ',num2str(years(ip2))];
    Xorig1 = [ones(size(TSMAX_4_10(ip1:ip2),1),1)  (TSMAX_4_10(ip1:ip2)) ]; %osservati
    [B_diff_per(:,ip)] = regress((yor(ip1:ip2)),Xorig1,alpha);
    bootdat = [(yor(ip1:ip2)),Xorig1];
    bootb = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
    bootCI(ip,:,:,:)=prctile(bootb,[2.5 97.5]);
end

%plot
y=[B1(2) B_detrend(2) B11(2) B21(2) B_smallest(2) B_largest(2) B_diff_per(2,:) ];
err=ones(1,length(y),2)*NaN;
err(1,:,1) = [B1(2)-bootCI1(1,2) B_detrend(2)-bootCI_detrend(1,2) B11(2)-bootCI11(1,2) B21(2)-bootCI21(1,2)  ...
    B_smallest(2)-bootCI_smallest(1,2) B_largest(2)-bootCI_largest(1,2) B_diff_per(2,:)-bootCI(:,1,2)'];

err(1,:,2) = [bootCI1(2,2)-B1(2) bootCI_detrend(2,2)-B_detrend(2) bootCI11(2,2)-B11(2) bootCI21(2,2)-B21(2)  ...
    bootCI_smallest(2,2)-B_smallest(2) bootCI_largest(2,2)-B_largest(2) bootCI(:,2,2)'-B_diff_per(2,:)];

figure(1); clf; 
hold on
for i = 1:length(y)
    h=bar(i,y(i));
    if (i>1)
        errorbar(i, y(:,i), err(:,i,1),err(:,i,2), 'LineStyle', 'none', ...
        'Color', [210 210 210]/255, 'LineWidth', 1);
    set(h,'FaceColor',[210 210 210]/255);
    errorbar(i, y(:,i), err(:,i,1),err(:,i,2), 'LineStyle', 'none', ...
        'Color', [0 0 0]/255, 'LineWidth', 1);
    else
        errorbar(i, y(:,i), err(:,i,1),err(:,i,2), 'LineStyle', 'none', ...
        'Color', [150 150 150]/255, 'LineWidth', 1);
    set(h,'FaceColor',[150 150 150]/255);
    errorbar(i, y(:,i), err(:,i,1),err(:,i,2), 'LineStyle', 'none', ...
        'Color', [0 0 0]/255, 'LineWidth', 1);
    end
end
xlabel('model calibration periods','Fontsize',font_size)
ylabel('TS_{max} coefficient','Fontsize',font_size)
gridxy([],y(:,1),'Color',[121 121 121]/255,'Linestyle','-');
set(gca,'XTick',1:size(xlabel_periods,1));
set(gca,'XTickLabel',xlabel_periods,'Fontsize',font_size);
htick=xticklabel_rotate([],45,[]);
set(gca,'FontSize',font_size)
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'stability_model_only_tsmax.eps']
print( gcf, '-depsc2', file ,'-painters')