CCGtl=FSpair(1).CCG.param.Timeline;

temp = [FSpair.CCG];
FSpairCCG_mega = [temp.Sessionwide];
FSpairCCGjm_mega = [temp.Sessionwide_JitterMean];
CCGsessionwide=FSpairCCG_mega-FSpairCCGjm_mega;
CCGdoub=conv2((CCGsessionwide+flip(CCGsessionwide,1)) / 2, 1/5*ones(5,1), 'same');

%%
fs=16;
pltopt=3;
switch pltopt
  case 1
    temppair=CCGdoub(CCGtl==0,:);
    xl=[-0.02 0.2]/5;
    yt=-0.03:0.01:0.03;
    ii=3;
    xlab='Excess Synchrony (Hz)';
    xtit='A';
  case 2
    temppair = [FSpair.SpikeCountCorr];
    xl=[-0.05 0.45];
    yt=-0.1:.1:0.5;
    ii=1;
    xlab='Spike Count Correlation';
    xtit='B';
  case 3
    temppair = [FSpair.CVCV_Spearman];
    xl=[-.05 .15];
    yt=-0.1:.05:0.2;
    ii=1;
    xlab='CV - CV  \rho_S_p_e_a_r_m_a_n';
    xtit='C';
end

figure('Position', [100 100 600 240])
xlabs=[];
tempcorr=[];
templab=[];
hold all
% plot([0 0], [0 7], '-', 'Color', [0 0 0 .1], 'LineWidth', 1.5)
for spi=1:6
  switch spi
    case 1
      p2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      ptit='sFS-sFS';
      pc=[.5 .5 .5];
    case 2
      p2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      ptit='grnsFS-grnsFS';
      pc=[0 .5 .5];
    case 3
      p2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      ptit='plnsFS-plnsFS';
      pc=[.5 .5 0];
    case 4
      p2p=find((strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS')) | ...
        (strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'sFS')));
      ptit='sFS-grnsFS';
      pc=.1+[.25 .5 .5];
    case 5
      p2p=find((strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS')) | ...
        (strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'sFS')));
      ptit='sFS-plnsFS';
      pc=.1+[.5 .5 .25];
    case 6
      p2p=find((strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS')) | ...
        (strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS')));
      ptit='grnsFS-plnsFS';
      pc=.1+[.25 .5 .25];
  end
  switch ii
    case 1
      p2p=p2p;
      stit='all pairs';
    case 2
      p2p = p2p(sametetFSpairs(p2p));
      stit='same tetrode pairs';
    case 3
      p2p = p2p(~sametetFSpairs(p2p));
      stit='not same tetrode pairs';
  end
  xlabs=[xlabs strcat(ptit, {' N='}, num2str(numel(p2p)))];
  tempcorr = [tempcorr temppair(p2p)];
  templab = [templab spi*ones(size(p2p))];
  
  plot(median(temppair(p2p)), spi, 's', 'MarkerEdgeColor', pc, 'MarkerFaceColor', pc, 'MarkerSize', 10)
  plot([prctile(temppair(p2p),25) prctile(temppair(p2p),75)], [spi spi], '-', 'Color', pc, 'LineWidth', 4)
  scatter(temppair(p2p), (spi+.4)*ones(size(p2p)), 8, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', pc, 'MarkerFaceAlpha', .5)
end
[p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(xl(2)-range(xl)/50*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', .5)
  plot(xl(2)+.4*range(xl)/50-range(xl)/50*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'MarkerSize', 5, 'LineWidth', .5)
end
ylim([0.5 6.5])
xlim(xl)
set(gca, 'XTick', yt, 'YTick', 1:6, 'YTickLabel', xlabs, 'YDir', 'reverse', 'FontSize', fs)
xlabel(xlab, 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS4', xtit))

%% CV rate correlations
temppre = [FSpair.CVrate_Spearman];
temppost = [FSpair.rateCV_Spearman];
xl=[-.15 .3];
yl=[0.5 9.5];
yt=-0.5:.1:0.5;
xlab='CV - rate  \rho_S_p_e_a_r_m_a_n';

figure('Position', [100 100 600 300])
xlabs=[];
tempcorr=[];
templab=[];
hold all
for spi=1:9
  switch spi
    case 1
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      ptit='sFS-sFS';
      pcol=[.5 .5 .5];
    case 2
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      ptit='sFS-grnsFS';
      pcol=[.5 .5 .5];
    case 3
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      ptit='sFS-plnsFS';
      pcol=[.5 .5 .5];
    case 4
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      ptit='grnsFS-sFS';
      pcol=[0 .5 .5];
    case 5
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      ptit='grnsFS-grnsFS';
      pcol=[0 .5 .5];
    case 6
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      ptit='grnsFS-plnsFS';
      pcol=[0 .5 .5];
    case 7
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      ptit='plnsFS-sFS';
      pcol=[.5 .5 0];
    case 8
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      ptit='plnsFS-grnsFS';
      pcol=[.5 .5 0];
    case 9
      pre2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      post2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      ptit='plnsFS-plnsFS';
      pcol=[.5 .5 0];
  end
  temppair = [temppre(pre2p) temppost(post2p)]';
  
  xlabs=[xlabs {strcat(ptit, ' N=', num2str(numel(pre2p)+numel(post2p)))}];
  tempcorr = [tempcorr; temppair];
  templab = [templab; spi*ones(numel(pre2p)+numel(post2p),1)];
  
  plot(median(temppair), spi, 's', 'MarkerEdgeColor', pcol, 'MarkerFaceColor', pcol, 'MarkerSize', 10)
  plot([prctile(temppair,25) prctile(temppair,75)], [spi spi], '-', 'Color', pcol, 'LineWidth', 4)
  scatter(temppair, (spi+.4)*ones(numel(pre2p)+numel(post2p),1), 8, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', pcol, 'MarkerFaceAlpha', .5)
end
xl=[-.15 .3];
yl=[0.5 9.5];
[p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(xl(2)-range(xl)/50*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', .5)
  plot(xl(2)+.4*range(xl)/50-range(xl)/50*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'MarkerSize', 5, 'LineWidth', .5)
end
set(gca, 'XTick', -.2:.1:.6, 'YTick', 1:spi, 'YTickLabel', xlabs, 'YDir', 'reverse', 'FontSize', fs)
xlim(xl)
ylim(yl)
xlabel(strcat('CV - Spike Count \rho_S_p_e_a_r_m_a_n'), 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS4D'))


%% FS subtypes waveform
PSUrepwf_mega = [FS.SpikeWaveform];
[PSUwfpeak, PSUwfpeakTind]=max(-PSUrepwf_mega,[],1);
PSUwftl=-PSUwfpeakTind+(1:40)';
PSUrepwfnorm_mega=PSUrepwf_mega./max(-PSUrepwf_mega,[],1);

yl=[-1 1.3];
figure('Position',[0 0 300 300])
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
  plot(1/30*PSUwftl(:,u2plot), PSUrepwfnorm_mega(:,u2plot), 'Color', [uc 0.2])
  % plot(1/30*PSUwftl(:,u2plot), -cumsum(PSUrepwfnorm_mega(:,u2plot),1), 'Color', [uc 0.2])
  text(1, yl(1)+.3*range(yl)-.1*range(yl)*spi, strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'Color', uc, ...
    'FontSize', 16, 'FontWeight', 'bold','HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
end
plot([-1/3 1], [0 0], 'k--', 'LineWidth', 0.5)
text(-1/3, 0, '0', 'FontSize', 16, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
xlim([-1/3 1])
ylim(yl)
axis('off')
print('-painters', '-dpdf', strcat(figPath, 'figS4Ei'))

%%
PSUwfpredipTind=NaN(1,numel(FS));
PSUwfpredip=NaN(1,numel(FS));
PSUwfahpTind=NaN(1,numel(FS));
PSUwfahp=NaN(1,numel(FS));
for psi=1:numel(FS)
  [PSUwfpredip(psi), PSUwfpredipTind(psi)] = max(PSUrepwf_mega(1:PSUwfpeakTind(psi),psi));
  
  [PSUwfahp(psi), tempTind] = max(PSUrepwf_mega(PSUwfpeakTind(psi):end,psi));
  PSUwfahpTind(psi) = PSUwfpeakTind(psi)-1 + tempTind;
end
if ~(isequal(PSUwfpredip, PSUrepwf_mega(PSUwfpredipTind + size(PSUrepwf_mega,1)*(0:numel(FS)-1))) ...
    && isequal(PSUwfahp, PSUrepwf_mega(PSUwfahpTind + size(PSUrepwf_mega,1)*(0:numel(FS)-1))))
  error('check algorithm')
end

PSUwfnormpredip = PSUrepwfnorm_mega(PSUwfpredipTind + size(PSUrepwfnorm_mega,1)*(0:numel(FS)-1));
PSUwfnormahp = PSUrepwfnorm_mega(PSUwfahpTind + size(PSUrepwfnorm_mega,1)*(0:numel(FS)-1));

PSUwfrpredipahp = PSUwfnormpredip./PSUwfnormahp;

PSUvalwf=1:numel(FS);


% rmpath(genpath('D:\Dropbox\Dropbox (Brown)\RESEARCH\analysis\chronux_2_11'));
[r,c]=find(PSUwftl==0);
[pks,locs]=findpeaks(reshape(PSUrepwf_mega,[],1));
tempTinds=repmat((1:40)',1,numel(FS));
valpks=pks(tempTinds(locs)>1 & tempTinds(locs)<40);
vallocs=locs(tempTinds(locs)>1 & tempTinds(locs)<40);
temppsinds=repmat(1:numel(FS),40,1);

troughinds=vallocs(PSUwftl(vallocs)>0);
troughpsi=temppsinds(troughinds);
[troughpsui, ipsi, ipsui]=unique(troughpsi,'first');
PSUwftroughTind=NaN(1,numel(FS));
PSUwftroughTind(troughpsui)=tempTinds(troughinds(ipsi));

PSUwfTpt=PSUwftroughTind-PSUwfpeakTind; % peak to trough time


figure('Position', [100 100 500 300])
% figure('Position', [100 100 750 300])
for whichwf=1:3
  switch whichwf
    case 1
      tempwf=1/30*PSUwfTpt;
      wftitle='Trough to Peak Time (ms) ';
      yl=[0.1 .3];
    case 2
      tempwf=PSUwfnormpredip;
      wftitle='Pre-Peak / Trough';
      yl=[-.1 0.6];
    case 3
      tempwf=PSUwfnormahp;
      wftitle='Peak / Trough';
      yl=[.2 .7];
    case 4
      tempwf=PSUwfrpredipahp;
      wftitle='Pre-Peak / Peak';
      yl=[-.4 1.3];
  end
  
  xlabs=[];
  tempwfs=[];
  templab=[];
  
  % figure('Position', [200*whichwf 200 200 300])
  subplot(1,3,whichwf)
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
    u2plot=PSUvalwf(ismember(PSUvalwf, u2plot));
    xlabs=[xlabs {utitle}];
    %   xlabs=[xlabs strcat(utitle, {' N='}, num2str(numel(u2plot)))];
    tempwfs = [tempwfs tempwf(u2plot)];
    templab = [templab spi*ones(1,numel(u2plot))];
    
    plot(spi, median(tempwf(u2plot)), 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 14)
    plot([spi spi], [prctile(tempwf(u2plot),25) prctile(tempwf(u2plot),75)], '-', 'Color', uc, 'LineWidth', 4)
    %   scatter((spi+.33)*ones(size(u2plot)), tempsta(u2plot), 8, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', uc, 'MarkerFaceAlpha', .5)
    
  end
  % yl=ylim;
  [p,tbl,stats] = kruskalwallis(tempwfs', templab', 'off');
  % text(.5, yl(1),sprintf('p=%.3f', p), 'FontSize', 16, 'HorizontalAlignment', 'left')
  % c=multcompare(stats,'CType','dunn-sidak','Display','off');
  c=multcompare(stats,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), yl(2)-range(yl)/34*[ci ci], 'k-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), yl(2)+range(yl)/34-range(yl)/34*[ci ci], 'k*', 'LineWidth', 1)
  end
  set(gca, 'XTick', 1:spi, 'XTickLabel', xlabs, 'XTickLabelRotation', 45, 'YTickLabelRotation', 90, 'FontSize', fs)
  xlim([0.5 spi+.5])
  ylim(yl)
  % ylim([yl(1)-.1*range(yl) yl(2)])
  ylabel(wftitle, 'FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'figS4Eii'))


%% FR stick IQR horizontal + histogram
tempvar=[FS.BaselineRate];
he=0:2:50;
xl=[0 25];
xlab='Spontaneous Rate (Hz)';

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [100 100 250 300])
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
  
  plot(median(tempvar(u2plot)), spi, 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 10)
  plot([prctile(tempvar(u2plot),25) prctile(tempvar(u2plot),75)], [spi spi], '-', 'Color', uc, 'LineWidth', 2.5)
  %   scatter((spi+.33)*ones(size(u2plot)), tempsta(u2plot), 8, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', uc, 'MarkerFaceAlpha', .5)
  
end
[p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(xl(2)-range(xl)/34*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', 1)
  plot(xl(2)+range(xl)/40-range(xl)/34*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'LineWidth', 1)
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
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(tempvar(u2plot), he, 'EdgeColor', 'none', 'FaceColor', uc, 'FaceAlpha', .5)
  end
end
set(gca, 'FontSize', fs, 'ActivePositionProperty', 'outerposition')
xlim(xl)
ylim([0 30])
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS4F'))

%% S2 PSD
ft=[8 15 30 60 120];
xl=[log10(ft(1)) log10(ft(end))];
yl=[8.25 13.25];

temp = [FS.SpikeSpectralPower];
temp = [temp.Pre1s];
temps2 = [temp.AllTrials];

figure('Position',[0 0 300 300])
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
  lw=2.5;
  shadedErrorBar(log10(FS(1).SpikeSpectralPower.param.frequency.Pre1s), mean(temps2(:,u2plot),2), ...
    std(temps2(:,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', uc, 'LineWidth', lw}, 1)
  
  text(xl(1), yl(2)-.1*range(yl)*(spi-1), strcat({' '}, utitle), ...
    'Color', uc, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
end
ylim(yl)
xlim(xl)
set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs) % , 'YTick', yl
xlabel('log Frequency (Hz)', 'FontSize', fs)
ylabel('Spike Spectral Power', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS4Gi'))


%% S2 stick IQR VERTICAL
figure('Position', [200 200 720 300])
for whichband=1:4
  switch whichband
    case 1
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=8 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<15;
      xlab='Power 8-15 Hz';
    case 2
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=15 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<30;
      xlab='Power 15-30 Hz';
    case 3
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=30 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<=55;
      xlab='Power 30-55 Hz';
    case 4
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=66 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<120;
      xlab='Power 66-120 Hz';
  end
  xl=[4 16];
  
  temps=geomean(temps2(sfcband,:), 1);
  
  
  xlabs=[];
  tempcorr=[];
  templab=[];
  
  subplot(1,4,whichband)
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
    tempcorr = [tempcorr temps(u2plot)];
    templab = [templab spi*ones(1,numel(u2plot))];
    
    plot(spi, median(temps(u2plot)), 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 14)
    plot([spi spi], [prctile(temps(u2plot),25) prctile(temps(u2plot),75)], '-', 'Color', uc, 'LineWidth', 4)
  end
  [p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
  disp(xlab); disp(p)
  c=multcompare(stats,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), xl(2)-range(xl)/34*[ci ci], 'k-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), xl(2)+range(xl)/34-range(xl)/34*[ci ci], 'k*', 'LineWidth', 1)
  end
  plot([0 7], [0 0], '-', 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
  set(gca, 'XTick', 1:spi, 'XTickLabel', xlabs, 'XTickLabelRotation', 45, 'FontSize', fs)%, 'YTick', 0:.2:2)
  xlim([0.5 spi+.5])
  ylim(xl)
  ylabel(xlab, 'FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'figS4Gii'))

