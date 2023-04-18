function [BA,ra,permrat,ratold,permratold]=ba_feedbackmodel2(BA,Z,Z2,totalarea,futureZ,futureP,option,capannual,nmem,firehistory,~)

Zveg=Z-Z2;
Zveg=(Zveg-mean(Zveg(1:51))./std(Zveg(1:51)));

% experimental function that builds in a negative feedback to burned area
nensemble=10000;
Zveg(6:end+5)=Zveg;
Z5=movmean(Zveg,[5 0]);
Z5=Z5(6:end);
Zveg=Zveg(6:end);


ZCMIPveg=futureZ-futureP;
ZCMIPveg=(ZCMIPveg-mean(ZCMIPveg(1:51))./std(ZCMIPveg(1:51)));

% experimental function that builds in a negative feedback to burned area
ZCMIPveg(6:end+5)=ZCMIPveg;
Z5CMIPveg=movmean(ZCMIPveg,[5 0]);
Z5CMIPveg=Z5CMIPveg(6:end);
ZCMIPveg=ZCMIPveg(6:end);


if nargin==7 % if you do not have fire history layer
    % SPINUP: assume some fraction of previous say 30-yrs have burnhistory
    % default it 0.1% of lands burn
    firehistory=0.001*ones(nmem,1);

    for i=2:length(BA)
        firehistory(1:nmem-1,i)=firehistory(2:nmem,i-1);
        firehistory(nmem,i)=BA(i-1)/totalarea;
    end
end


% reduce available area to burn as a function of fire history
addrat=0;
for i=1:length(BA)
    [rat,addrat,firehistory]=firehistorybufferadv(totalarea,firehistory,option,nmem,Zveg(i),Z5(i),addrat,0.1);
    firehistory(end+1)=BA(i)/totalarea;
    firehistory=firehistory(2:end);
    rrr(i)=rat+addrat;
    permoff(i)=addrat;
end
rat=rrr;
%Z2=Z;

k2=table(Z',Z2',log10(BA./(1-rat)'));
% statistics about the model fit
o=fitlm(k2);
mm=predict(o,k2);mm=(10.^mm)./(1-rat)';

%r,p,BA,ra,permrat,ratold,permratold
% residuals and constrained residuals
residual_log=(log10(BA)-log10(mm));
sigma=std(residual_log);

% save permoff firehistory and rrr
permratold=permoff;
ratold=rat;

reset1950=0;
if reset1950
    % if you want to reset burned area starting at 1950
    firehistory=firehistory1950;
    clear BA rat
    startyr=1;
else
    startyr=72; % this should be years since 1950 - basically the starting year
end
ra(startyr-length(rat):startyr-1,1)=rat;
BA(startyr-length(rat):startyr-1,1)=BA;
if option==0 startyr=1;end  % this is for no-feedback model

% starting in 1950 iterative forward
%Z5(6:156)=movmean(futureZ(1:151),[5 0]);
addrat=0;

%k3=table(futureZ',futureP');
% statistics about the model fit
%o=fitlm(k2);
k3=[futureZ futureP];

mm=predict(o,k2);mm=(10.^mm)./(1-rat)';

initfhist=firehistory;

%ki=table(futureZ,futureP);


for kk=1:nensemble
    if mod(kk,100) == 0
        fprintf('At iteration %d...\n',kk);
    end
    firehistory=initfhist;addrat=0;
    for i=startyr:151
        [rat,addrat,firehistory]=firehistorybufferadv(totalarea,firehistory,option,nmem,ZCMIPveg(i),Z5CMIPveg(i),addrat,.1);
        noise=normrnd(0,sigma,1,1);
        %ki=table(futureZ(i),futureP(i));
        
        BA(i,kk)=10.^(predict(o,k3(i,:))+noise).*(1-rat);

        % limit such that the burned area annually can not
        % exceed more than X% of existing land
        BA(i,kk)=min(BA(i,kk),(1-rat)*totalarea*capannual);
        ra(i,kk)=rat;
        permrat(i,kk)=addrat;
        firehistory(end+1)=BA(i,kk)/totalarea;
        firehistory=firehistory(2:end);

    end
end