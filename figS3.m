%% pre 1s all trials vs sessionwide comparison

tpre = FS(1).ISI.param.Sliding250msTimeline==-125;

temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
isipeakvec=[temp.ISIpeak];

temp = [FS.ISI];
temp=[temp.Pre100ms_2ndST];
temp=[temp.AllTrials];
isipeakvec2ndst=[temp.ISIpeak];

temp = [FS.ISI];
temp=[temp.Sessionwide];
isipeakvecsw=[temp.ISIpeak];

temp = [FS.ISI];
temp=[temp.Sliding250ms];
temp=[temp.AllTrials];
isipeakvecsliding=[temp.ISIpeak];


temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
cvvec=[temp.CV];

temp = [FS.ISI];
temp=[temp.Pre100ms_2ndST];
temp=[temp.AllTrials];
cvvec2ndst=[temp.CV];

temp = [FS.ISI];
temp=[temp.Sessionwide];
cvvecsw=[temp.CV];

temp = [FS.ISI];
temp=[temp.Sliding250ms];
temp=[temp.AllTrials];
cvvecsliding=[temp.CV];

fs=8;
ft=[8 15 30 60 120];
figure('Position', [100 100 300 450])
for whichparam=1:6
  switch whichparam
    case 1
      xtemp=isipeakvec;
      xlab='Pre 1s';
      ytemp=isipeakvecsw;
      ylab='Sessionwide';
      titlab='ISI Peak (ms)';
      yl=[0 50];
      yt=yl;
    case 3
      xtemp=isipeakvecsliding(tpre,:);
      xlab='Pre 250ms';
      ytemp=isipeakvec;
      ylab='Pre 1s';
      titlab='ISI Peak (ms)';
      yl=[0 50];
      yt=yl;
    case 5
      xtemp=isipeakvec2ndst;
      xlab='Pre 100ms 2^n^d ST';
      ytemp=isipeakvecsw;
      ylab='Sessionwide';
      titlab='ISI Peak (ms)';
      yl=[0 50];
      yt=yl;
    case 2
      xtemp=cvvec;
      xlab='Pre 1s';
      ytemp=cvvecsw;
      ylab='Sessionwide';
      titlab='CV';
      yl=[.8 2];
      yt=0:.5:5;
    case 4
      xtemp=cvvecsliding(tpre,:);
      xlab='Pre 250ms';
      ytemp=cvvec;
      ylab='Pre 1s';
      titlab='CV';
      yl=[.6 1.6];
      yt=0:.5:5;
    case 6
      xtemp=cvvec2ndst;
      xlab='Pre 100ms 2^n^d ST';
      ytemp=cvvecsw;
      ylab='Sessionwide';
      titlab='CV';
      yl=[.8 2];
      yt=0:.5:5;
  end
  
  subplot(3,2,whichparam)
  hold all
  for spi=1:3
    switch spi
      case 1
        u2plot=sFS;
        utitle='sFS';
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
        uc=[0 .5 .5];
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
        uc=[.5 .5 0];
    end
    scatter(xtemp(u2plot), ytemp(u2plot), 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', uc, 'MarkerFaceAlpha', 0.4)
    if whichparam==1
      text(xl(1), yl(2)-.1*range(yl)*(spi-1), strcat({' '},utitle), ...
        'Color', uc, 'Fontsize', fs*7/8, 'Fontweight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    end
  end
  text(yl(2),yl(1), strcat('\rho_S_p_e_a_r_m_a_n=',sprintf('%.3f', corr(xtemp', ytemp', 'Type', 'Spearman'))), ...
    'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', fs*7/8)
  plot(yl, yl, 'Color',[0 0 0 0.1], 'LineWidth', 1.5)
  axis([yl yl])
  set(gca, 'FontSize', fs);
  xlabel(xlab, 'FontSize', fs);
  ylabel(ylab, 'FontSize', fs);
  if whichparam<=2
    title(titlab, 'FontSize', fs)
  else
    title(' ', 'FontSize', fs)
  end
end
print('-painters', '-dpdf', strcat(figPath, 'figS3AB'))

%% CV and rate relationship
fs=16;

sponratevec = [FS.BaselineRate];

temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
cvvec=[temp.CV];

xl=[0 25];
yl=[0.8 1.6];

figure('Position', [0 0 300 300])
hold all
for spi=0:3
  switch spi
    case 0
      u2plot=supFS;
      utitle='supFS';
      uc=[0 0 0];
    case 1
      u2plot=sFS;
      utitle='sFS';
      uc=[.5 .5 .5];
    case 2
      u2plot=grnsFS;
      utitle='grnsFS';
      uc=[0 .5 .5];
    case 3
      u2plot=plnsFS;
      utitle='plnsFS';
      uc=[.5 .5 0];
  end
  scatter(sponratevec(u2plot), cvvec(u2plot), 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', uc, 'MarkerFaceAlpha', 0.4)
  if spi>0
    text(xl(2), yl(2)-.1*range(yl)*(spi-1), strcat({' '},utitle, ' (N=', num2str(numel(u2plot)), ')'), ...
      'Color', uc, 'Fontsize', fs*7/8, 'Fontweight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
  end
end
utmquintile=prctile(sponratevec,0:20:100);
utmquintmean=zeros(5,1);
utmquinterr=zeros(5,1);
cvsqquintmean=zeros(5,1);
cvsqquinterr=zeros(5,1);
for ii=1:5
  quintu= (sponratevec>=utmquintile(ii) &  sponratevec<=utmquintile(ii+1));
  utmquintmean(ii)=mean(sponratevec(quintu));
  utmquinterr(ii)=std(sponratevec(quintu))/sqrt(nnz(quintu));
  cvsqquintmean(ii)=mean(cvvec(quintu));
  cvsqquinterr(ii)=std(cvvec(quintu))/sqrt(nnz(quintu));
end
errorbar(utmquintmean, cvsqquintmean, cvsqquinterr, cvsqquinterr, utmquinterr, utmquinterr, 'k.-', 'LineWidth', 4, 'CapSize', 0)
text(xl(1),yl(1), strcat('\rho_S_p_e_a_r_m_a_n=', sprintf('%.3f', corr(sponratevec', cvvec', 'Type', 'Spearman'))), ...
  'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 14)
set(gca, 'FontSize', fs)
ylim(yl)
xlim(xl)
xlabel('Rate (Hz)', 'FontSize', fs)
ylabel('CV', 'FontSize', fs)
title(' ', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS3C'))

%% CV stick IQR horizontal + histogram
fs=8;
temp = [FS.SpearmanCorr];
tempvar = [temp.CVrate];

he=-1:.08:1;
xl=[-0.15 .5];
xlab='CV-Rate \rho_S_p_e_a_r_m_a_n';
  
xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [0 200 125 150])
subplot(5,1,1)
hold all
for spi=1:3
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
      uc=[.5 .5 .5];
    case 2
      u2plot=grnsFS;
      utitle='grnsFS';
      uc=[0 .5 .5];
    case 3
      u2plot=plnsFS;
      utitle='plnsFS';
      uc=[0.5 0.5 0];
  end
  
  xlabs=[xlabs {utitle}];
%   xlabs=[xlabs strcat(utitle, {' N='}, num2str(numel(u2plot)))];
  tempcorr = [tempcorr tempvar(u2plot)];
  templab = [templab spi*ones(1,numel(u2plot))];
  
  plot(median(tempvar(u2plot)), spi, 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 5)
  plot([prctile(tempvar(u2plot),25) prctile(tempvar(u2plot),75)], [spi spi], '-', 'Color', uc, 'LineWidth', 1.25)
%   scatter((spi+.33)*ones(size(u2plot)), tempsta(u2plot), 8, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', uc, 'MarkerFaceAlpha', .5)
  
end
[p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
disp(p)
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(xl(2)-range(xl)*0.05*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', .25)
  plot(xl(2)+range(xl)*0.02-range(xl)*0.05*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'MarkerSize', 5, 'LineWidth', .25)
end
axis off
set(gca, 'XTick', [], 'YTick', [], 'FontSize', fs)%, 'ActivePositionProperty', 'outerposition') 
ylim([0.5 spi+.5])
xlim(xl)

subplot(5,1,2:5)
hold all
for spi=1:3
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
      uc=[.5 .5 .5];
    case 2
      u2plot=grnsFS;
      utitle='grnsFS';
      uc=[0 .5 .5];
    case 3
      u2plot=plnsFS;
      utitle='plnsFS';
      uc=[0.5 0.5 0];
  end
  
  if spi==1  
    hc=histcounts(tempvar(u2plot), he);
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',2)
  else
    histogram(tempvar(u2plot), he, 'EdgeColor', 'none', 'FaceColor', uc, 'FaceAlpha', .5)
  end
end
set(gca, 'FontSize', fs, 'ActivePositionProperty', 'outerposition')
xlim(xl)
ylim([0 35])
xlabel(xlab, 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS3D'))

%% schematic: inverse ISI and spike spectral power relationship
fs=16;
addpath(genpath('D:\Dropbox\Dropbox (Brown)\RESEARCH\analysis\chronux_2_11'));

tapers=[3 5]; 
Fs=1000;
pad=0; 
err=0;
trialave=0;
fpass=[0 500];

dt=1/Fs;
t=0:dt:.25; % time grid for prolates
N=length(t); % number of points in grid for dpss
nfft=max(2^(nextpow2(N)+pad),N); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass);
tapers=dpsschk(tapers,N,Fs); % check tapers


% isivec=.005:.005:.245;
isivec=0:dt:.25;
% isivec=10.^-(log10(1/0.25):0.01:log10(1/0.001));
S2mat=NaN(length(isivec), length(f));

for ii=1:length(isivec)
  tempisi=isivec(ii);
  
  [J2,Msp2,Nsp2]=mtfftpt(.125+.5*[-tempisi; tempisi],tapers,nfft,t,f,findx);
  S2=squeeze(mean(conj(J2).*J2,2)); % spectrum data 2
  S2mat(ii,:)=S2;
end

rmpath(genpath('D:\Dropbox\Dropbox (Brown)\RESEARCH\analysis\chronux_2_11'));

ft=[8 15 30 60 120];
figure('Position', [0 0 400 300])
hold all
h=pcolor(log10(1./isivec),log10(f),log10(S2mat)');
set(h, 'EdgeColor', 'none')
plot(log10([f(2) f(end)]), log10(1/20*[f(2) f(end)]), 'r-')
for ii=1:4
  plot(log10([f(2) f(end)]), log10(1/ii*[f(2) f(end)]), 'w-')
  plot(log10([f(2) f(end)]), log10(ii*[f(2) f(end)]), 'w-')
end
set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', log10(ft), 'YTickLabel', ft, 'FontSize', fs)
axis(log10([ft(1) ft(end) ft(1) ft(end)]))
% axis(log10([8 80 8 80]))
xlabel('log ISI^-^1 (s^-^1)', 'FontSize', fs)
ylabel('log Frequency (Hz)', 'FontSize', fs)
caxis([0 1.5])
colorbar
print('-painters', '-dpdf', strcat(figPath, 'figS3E'))

%%
fs=8;
ft=[8 15 30 60 120];

temprho = [FS.SpearmanCorr];

figure('Position', [100 100 450 750])
for pltopt=1:5  
  switch pltopt
    case 1
      tempxvec=FS(1).SpearmanCorr.frequency;
      tempcorr=[temprho.LFPSpikeSpectralPower];
      temppcorr=[temprho.pLFPSpikeSpectralPower];
      plttitle='LFP-Spike Spectral Power        ';
      xlab='Frequency (Hz)';
      yl=[-.15 .15];
    case 3
      tempxvec=FS(1).SpearmanCorr.frequency;
      tempcorr=[temprho.rateSpikeSpectralPower];
      temppcorr=[temprho.prateSpikeSpectralPower];
      plttitle='Rate-Spike Spectral Power    ';
      xlab='Frequency (Hz)';
      yl=[.5 .8];
    case 2
      tempxvec=1000./(FS(1).ISI.param.NormProbtimeline);
      tempcorr=[temprho.rateISIrate];
      temppcorr=[temprho.prateISIrate];
      plttitle='Rate-ISI Rate';
      xlab='ISI^-^1 (s^-^1)';
      yl=[0 .6];
    case 5
      tempxvec=FS(1).SpearmanCorr.frequency;
      tempcorr=[temprho.CVSpikeSpectralPower];
      temppcorr=[temprho.pCVSpikeSpectralPower];
      plttitle='CV-Spike Spectral Power    ';
      xlab='Frequency (Hz)';
      yl=[-.15 .35];
    case 4
      tempxvec=1000./(FS(1).ISI.param.NormProbtimeline);
      tempcorr=[temprho.CVISIrate];
      temppcorr=[temprho.pCVISIrate];
      plttitle='CV-ISI Rate';
      xlab='ISI^-^1 (s^-^1)';
      yl=[-.15 .35];
  end
  
  xl=log10([ft(1) ft(end)]);
  subplot(5,3,3*(pltopt-1)+1)
  hold all
  for spi=1:3
    switch spi
      case 1
        u2plot=sFS;
        utitle='sFS';
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
        uc=[0 .5 .5];
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
        uc=[0.5 0.5 0];
    end
    shadedErrorBar(log10(tempxvec), nanmean(tempcorr(:,u2plot),2), ...
      nanstd(tempcorr(:,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', uc, 'LineWidth', 1.25},1)
    
    if pltopt==1
      text(xl(1), yl(1)+.2*range(yl)-.1*range(yl)*(spi-1), strcat({' '}, utitle), ...
        'Color', uc, 'FontWeight', 'bold', 'FontSize', 7, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    end
  end
  % yl=ylim;
  plot(log10([ft(1) ft(end)]), [0 0], 'Color', [0 0 0 0.1], 'LineWidth', .75)
  ylim(yl)
  xlim(xl)
  if contains(plttitle, 'CV')
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YDir', 'reverse', 'FontSize', fs)
  else
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs)
  end
  xlabel(strcat('log', {' '}, xlab), 'FontSize', fs)
  ylabel('\rho_S_p_e_a_r_m_a_n', 'FontSize', fs)
  title(plttitle, 'FontSize', fs)%, 'FontWeight', 'normal')
  
  
  sigposcorr=tempcorr>0 & temppcorr<0.05;
  subplot(5,3,3*(pltopt-1)+2)
  hold all
  plot(log10([tempxvec(1) tempxvec(end)]), [2.5 2.5], 'Color', [0 0 0 0.1], 'LineWidth', .75)
  for spi=1:3
    switch spi
      case 1
        u2plot=sFS;
        utitle='sFS';
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
        uc=[0 .5 .5];
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
        uc=[0.5 0.5 0];
    end
    plot(log10(tempxvec), 100*mean(sigposcorr(:,u2plot),2), 'Color', uc, 'LineWidth', 1.25)
  end
  set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs) %, 'XGrid', 'on', 'YGrid', 'on')
  xlabel(strcat('log', {' '}, xlab), 'FontSize', fs)
  ylabel('% Significant Units', 'FontSize', fs)
  xlim(log10([ft(1) ft(end)]))
  ylim([0 100])
  title('Positive Correlation', 'FontSize', fs)
  
  signegcorr=tempcorr<0 & temppcorr<0.05;
  subplot(5,3,3*pltopt)
  hold all
  plot(log10([tempxvec(1) tempxvec(end)]), [2.5 2.5], 'Color', [0 0 0 0.1], 'LineWidth', .75)
  for spi=1:3
    switch spi
      case 1
        u2plot=sFS;
        utitle='sFS';
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
        uc=[0 .5 .5];
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
        uc=[0.5 0.5 0];
    end
    plot(log10(tempxvec), 100*mean(signegcorr(:,u2plot),2), 'Color', uc, 'LineWidth', 1.25)
  end
  set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs) %, 'XGrid', 'on', 'YGrid', 'on')
  xlabel(strcat('log', {' '}, xlab), 'FontSize', fs)
  ylabel('% Significant Units', 'FontSize', fs)
  xlim(log10([ft(1) ft(end)]))
  ylim([0 100])
  title('Negative Correlation', 'FontSize', fs)  
end
print('-painters', '-dpdf', strcat(figPath, 'figS3F'))
