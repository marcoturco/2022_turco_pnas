function dataout=scale_base_period(datain,base_period,years_series)

[~,~,iok]=intersect(base_period, years_series);
avg=nanmean(datain(iok),1);
devstd=std(datain(iok),'omitnan');
dataout=(datain-avg)/devstd;
%dataout=datain-avg;

    




     