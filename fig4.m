load(strcat(dataPath, 'FSpair.mat'))
sametetFSpairs = [FSpair.issametet];

%% CCG all pairs

CCGtl=FSpair(1).CCG.param.Timeline;

temp = [FSpair.CCG];
FSpairCCG_mega = [temp.Sessionwide];
FSpairCCGjm_mega = [temp.Sessionwide_JitterMean];
CCGsessionwide=FSpairCCG_mega-FSpairCCGjm_mega;
CCGdoub=conv2((CCGsessionwide+flip(CCGsessionwide,1)) / 2, 1/5*ones(5,1), 'same');

yl=[-.005 .018];
yt=-0.03:0.01:0.03;
xl=[-95 95];

figure('Position',[0 0 900 300])
for spi=1:3
  switch spi
    case 1
      p2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      ptit='sFS Pairs';
      pc=[.5 .5 .5];
    case 2
      p2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      ptit='grnsFS Pairs';
      pc=[0 .5 .5];
    case 3
      p2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      ptit='plnsFS Pairs';
      pc=[.5 .5 0];
  end
  p2p = p2p(~sametetFSpairs(p2p));
  stit='not same tetrode pairs';
  subplot(1,3,spi)
  hold all
  shadedErrorBar(CCGtl, mean(CCGdoub(:,p2p), 2), std(CCGdoub(:,p2p),0,2)/sqrt(2*numel(p2p)), {'Color', pc, 'LineWidth', 2.5}, 1)
  plot([-100 100], [0 0], '-', 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
  xlim(xl)
  ylim(yl)
  title(strcat(ptit, ' (N=', num2str(numel(p2p)), ')'), 'FontSize', fs)
  if spi==2
    xlabel('Time Lag from Spike (ms)', 'FontSize', fs)
  end
  if spi==1
    ylabel('Excess Rate (Hz)', 'FontSize', fs)
  end
  set(gca, 'XTick', -100:50:100, 'FontSize', fs, 'YTick', yt, 'YTickLabel', sprintf('%g\n', yt))
end
print('-painters', '-dpdf', strcat(figPath, 'fig4A'))

%% excess synch count VERTICAL
fs=16;

tempsynch=CCGdoub(CCGtl==0,:);
xl=[-0.004 0.025];
yt=-0.03:0.01:0.03;
ii=3;
xlab='Excess Synchrony (Hz)';

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [0 0 200 300])
hold all
for spi=1:3
  switch spi
    case 1
      p2p=find(strcmp({FSpair.subtype_FS1}, 'sFS') & strcmp({FSpair.subtype_FS2}, 'sFS'));
      ptit='sFS';
      pc=[.5 .5 .5];
    case 2
      p2p=find(strcmp({FSpair.subtype_FS1}, 'grnsFS') & strcmp({FSpair.subtype_FS2}, 'grnsFS'));
      ptit='grnsFS';
      pc=[0 .5 .5];
    case 3
      p2p=find(strcmp({FSpair.subtype_FS1}, 'plnsFS') & strcmp({FSpair.subtype_FS2}, 'plnsFS'));
      ptit='plnsFS';
      pc=[.5 .5 0];
  end
  p2p = p2p(~sametetFSpairs(p2p));
  stit='not same tetrode pairs';
  
  xlabs=[xlabs {ptit}];
  tempcorr = [tempcorr tempsynch(p2p)];
  templab = [templab spi*ones(size(p2p))];
  
  plot(spi, median(tempsynch(p2p)), 's', 'MarkerEdgeColor', pc, 'MarkerFaceColor', pc, 'MarkerSize', 14)
  plot([spi spi], [prctile(tempsynch(p2p),25) prctile(tempsynch(p2p),75)], '-', 'Color', pc, 'LineWidth', 4)
  scatter((spi+.4)*ones(size(p2p)), tempsynch(p2p), 14, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', pc, 'MarkerFaceAlpha', .5)
  
end
[p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(c(sigp(ci),1:2), xl(2)-range(xl)/34*[ci ci], 'k-', 'LineWidth', 1)
  plot(mean(c(sigp(ci),1:2)), xl(2)+range(xl)/34-range(xl)/34*[ci ci], 'k*', 'LineWidth', 1)
end
plot([0 7], [0 0], '-', 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
set(gca, 'XTick', 1:spi, 'XTickLabel', xlabs, 'XTickLabelRotation', 45, 'FontSize', fs, 'YTick', yt, 'YTickLabel', sprintf('%g\n', yt))
xlim([0.5 spi+.5])
ylim(xl)
ylabel(xlab, 'FontSize', fs)
% title('CCG Peak', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig4B'))


%% STA
temp = [FS.SpikeTriggeredAverageLFP];
temp = [temp.Pre1s];
STAmat = [temp.AllTrials];

xl=[-100 100];
yl=[-7 7];
figure('Position', [300 300 300 300])
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
  plot(FS(1).SpikeTriggeredAverageLFP.param.Timeline, mean(STAmat(:,u2plot),2), 'Color', uc, 'LineWidth',lw)
  text(xl(1), yl(1)+.31*range(yl)-.1*range(yl)*spi, strcat({' '}, utitle), 'Color', uc, 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom')
end
set(gca, 'FontSize', fs) %, 'ActivePositionProperty', 'outerposition')
xlim(xl)
ylim(yl)
% ylim([-10 10])
title('Spike Triggered Average', 'FontSize',fs)
ylabel('LFP (\muV)', 'FontSize',fs)
xlabel('Time from Spike (ms)', 'FontSize',fs)

print('-painters', '-dpdf', strcat(figPath, 'fig4Ci'))


%% STA stick IQR horizontal + histogram
STA0ind= find(FS(1).SpikeTriggeredAverageLFP.param.Timeline==0);

tempsta=1/5*(mean(STAmat(STA0ind-2:STA0ind+2,:),1) - mean(STAmat(STA0ind-7:STA0ind-3,:),1));
xl=1/5*[-6 3];
xlab='\DeltaLFP at T=0ms (\muV/ms)';
he=1/5*[-50:.5:50];

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [300 200 250 300])
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
  tempcorr = [tempcorr tempsta(u2plot)];
  templab = [templab spi*ones(1,numel(u2plot))];
  
  plot(median(tempsta(u2plot)), spi, 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 10)
  plot([prctile(tempsta(u2plot),25) prctile(tempsta(u2plot),75)], [spi spi], '-', 'Color', uc, 'LineWidth', 2.5)
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
    hc=histcounts(tempsta(u2plot), he);
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(tempsta(u2plot), he, 'EdgeColor', 'none', 'FaceColor', uc, 'FaceAlpha', .5)
  end
end
set(gca, 'FontSize', fs, 'ActivePositionProperty', 'outerposition')
xlim(xl)
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig4Cii'))

%% SFC PSD
temp = [FS.SpikeFieldCoherence];
temp = [temp.Pre1s];
SFCmat = [temp.AllTrials];

ft=[8 15 30 60 120];
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
  plot(log10(FS(1).SpikeSpectralPower.param.frequency.Pre1s), mean(log10(SFCmat(:,u2plot)),2), 'Color', uc, 'LineWidth', lw)
  text(log10(ft(end)), yl(2)-.1*range(yl)*(spi-1), strcat(utitle), 'Color', uc, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
end
yl=[-2.1 -1.2];
ylim(yl)
xlim([log10(ft(1)) log10(ft(end))])
set(gca, 'YTick', yl, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs) %, 'ActivePositionProperty', 'outerposition')
xlabel('log Frequency (Hz)', 'FontSize', fs)
ylabel('log Power', 'FontSize', fs)
title('Spike Field Coherence', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig4Di'))

%% SFC stick IQR VERTICAL
figure('Position', [200 200 720 300])
for whichband=1:4
  switch whichband
    case 1
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=8 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<15;
      xlab='log Power 8-15 Hz';
      sfctit='Alpha';
      xl=[-1.8 -1.2];
    case 2
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=15 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<30;
      xlab='log Power 15-30 Hz';
      sfctit='Beta';
      xl=[-2.05 -1.4];
    case 3
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=30 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<=55;
      xlab='log Power 30-55 Hz';
      sfctit='Gamma';
      xl=[-2.05 -1.6];
    case 4
      sfcband=FS(1).SpikeSpectralPower.param.frequency.Pre1s>=66 & FS(1).SpikeSpectralPower.param.frequency.Pre1s<120;
      xlab='log Power 66-120 Hz';
      sfctit='High-Gamma';
      xl=[-2.05 -1.75];
  end
  
  tempsfc=log10(geomean(SFCmat(sfcband,:), 1));
  
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
    tempcorr = [tempcorr tempsfc(u2plot)];
    templab = [templab spi*ones(1,numel(u2plot))];
    
    plot(spi, median(tempsfc(u2plot)), 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 14)
    plot([spi spi], [prctile(tempsfc(u2plot),25) prctile(tempsfc(u2plot),75)], '-', 'Color', uc, 'LineWidth', 4)
    
  end
  [p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
  c=multcompare(stats,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), xl(2)-range(xl)/34*[ci ci], 'k-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), xl(2)+range(xl)/34-range(xl)/34*[ci ci], 'k*', 'LineWidth', 1)
  end
  set(gca, 'XTick', 1:spi, 'XTickLabel', xlabs, 'XTickLabelRotation', 45, 'YTickLabelRotation', 90, 'FontSize', fs, 'YTick', -2:.2:-1)
  xlim([0.5 spi+.5])
  ylim(xl)
  ylabel(xlab, 'FontSize', fs)
  title(sfctit, 'FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'fig4Dii'))
