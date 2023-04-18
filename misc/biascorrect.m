function [cmip6]=biascorrect(Z,cmip6,nmod);

[~,~,iok_gcm]=intersect(1995:2014, 1950:2100);

% loop over models
for model=1:nmod
    [Z1(:,model),ZZ1(:,model)]=bias_correction_ecdf...
        (Z(end-29:end),cmip6(42:71,model),cmip6(:,model));
%    [Z1(:,model),ZZ1(42:71,model)]=bias_correction_ecdf...
%        (Z(end-29:end),cmip6(42:71,model),cmip6(42:71,model));
%    [Z1(:,model),ZZ1(72:101,model)]=bias_correction_ecdf...
%        (Z(end-29:end),cmip6(42:71,model),cmip6(72:101,model));
%    [Z1(:,model),ZZ1(102:151,model)]=bias_correction_ecdf...
%        (Z(end-29:end),cmip6(42:71,model),cmip6(102:151,model));
 %   avg=nanmean(ZZ1(iok_gcm,model),1);
 %   devstd=std(ZZ1(iok_gcm,model),'omitnan');
 %   ZZ1(1:end,model)=(ZZ1(1:end,model)-avg)/devstd;
end


cmip6=ZZ1;
