%%
phasecol=[
  0.5 0 1
  0 0 1
  0 0.7 0
  1 0.5 0
  1 0 0
  ];

xr=0;
replacenan=1;

temp=[FS.SpikeTiming];
temp=cat(1,temp.DeltaHitRate);

figure('Position', [0 100 1200 300])
for spi=0:3
  yl=[-0.12 0.12];
  
  switch spi
    case 0
      u2plot=1:numel(FS);
      utitle='FS';
      uc=[.5 .5 .5];
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
      uc=[0.5 0.5 0 .5];
  end
  
  fslastprespkthr=temp(u2plot,:)';
  
  subplot(1,4,1+spi)
  hold all
  plot([-25 0], [0 0], 'k--')
  plot(FS(1).SpikeTiming.DeltaHitRateParam, nanmean(fslastprespkthr,2), 'ko-', 'LineWidth', 2.5)
  for fi=1:size(fslastprespkthr,1)
    errorbar(-27.5+5*fi, nanmean(fslastprespkthr(fi,:)), nanstd(fslastprespkthr(fi,:))/sqrt(numel(u2plot)), ...
      'o-','Color', phasecol(fi,:), 'MarkerFaceColor', phasecol(fi,:), 'LineWidth', 4)
  end
  
  if replacenan
    fslastprespkthr(isnan(fslastprespkthr))=0;
    disp('replace nan')
  else % remove nan rows
    u2plot=u2plot(~any(isnan(fslastprespkthr), 1));
    fslastprespkthr(:,any(isnan(fslastprespkthr),1))=[];
    disp('remove nan rows')
  end
  
  tempmat=fslastprespkthr';
  [pfried, tblfried, statsfried]=friedman(tempmat, 1, 'off');
  
  c=multcompare(statsfried,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(-27.5+5*c(sigp(ci),1:2), .1-.05*range(yl)*[ci ci], 'k-', 'LineWidth', 1)
    plot(-27.5+5*mean(c(sigp(ci),1:2)), .1+.02*range(yl)-.05*range(yl)*[ci ci], 'k*', 'LineWidth', 1)
  end
  
  phasepvals=NaN(size(fslastprespkthr,1),1);
  for fi=1:size(fslastprespkthr,1)
    phasepvals(fi)=signrank(fslastprespkthr(fi,:));
  end
  phaseh=phasepvals<0.05/5;
  if nnz(phaseh)>0 && pfried<0.05
    plot(-27.5+5*find(phaseh==1), 0.08, 'r*', 'LineWidth', 1)
  end
  
  if xr==1
    text(-25, .11, sprintf(' p=%0.3f', pfried), 'FontSize', 16, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
    set(gca, 'XTick', -25:5:0, 'FontSize', fs, 'XDir', 'reverse')%, 'XTickLabelRotation', 45
  else
    text(0, .11, sprintf(' p=%0.3f', pfried), 'FontSize', 16, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
    set(gca, 'XTick', -25:5:0, 'FontSize', fs)%, 'XDir', 'reverse', 'XTickLabelRotation', 45
  end
  
  xlim([-25 0])
  ylim(yl)
  
  if spi==0
    ylabel('\DeltaHit Rate', 'FontSize', fs)
  end
  xlabel(' ', 'FontSize', fs)
  annotation('textbox', [0.25 0 0.5 0.1], 'String', 'Last Pre Spike Timing (ms)', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'normal', 'FontSize', fs, 'EdgeColor', 'none')
  %   title(' ', 'FontSize', fs)
  title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'fig7B'))

%% plot sPSTH5DP
temp = [FS.SpikeTiming];
temp = [temp.PeriDP];
temp=[temp.Sliding5ms]';

figure('Position', [0 100 1200 300])
for spi=0:3
  switch spi
    case 0
      u2plot=1:numel(FS);
      utitle='FS';
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
  subplot(1,4,spi+1)
  hold all
  plot([-250 250], [0.5 0.5], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
  %   plot([-12.5 -12.5], [0 1], 'Color', [0 1 0 .5], 'LineWidth', 1.5)
  for fi=1:size(phasecol,1)
    plot(-27.5+5*[fi fi], [0 1], 'Color', [phasecol(fi,:) .5])
  end
  
  shadedErrorBar(FS(1).SpikeTiming.PeriDP.param.Timeline, mean(temp(u2plot,:),1), ...
    std(temp(u2plot,:),0,1)/sqrt(numel(u2plot)), {'k-', 'LineWidth',2.5}, 1)
  
  xlim([-50 50])
  ylim([.48 .52])
    
  set(gca, 'XTick', -250:25:250, 'YTick', .48:.02:.52, 'FontSize', fs)
  % if spi==2
  % xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
  % end
  if spi==0
    ylabel('Detect Probability', 'FontSize', fs)
  end
  
  xlabel(' ', 'FontSize', fs)
  annotation('textbox', [0.25 0 0.5 0.1], 'String', 'Time from Sensory Onset (ms)', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'normal', 'FontSize', fs, 'EdgeColor', 'none')
  
  title(' ', 'FontSize', fs)
  %   title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'fig7C'))

%%
ft=[8 15 30 60 120];
xt=log10(ft);
xl=log10([ft(1) ft(end)]);
yl=[5*10^-7 6*10^-6];
xlab='log Frequency (Hz)';
ylab='Power';

f_DP5=FS(1).SpikeTiming.PeriDP.param.frequency;
temp = [FS.SpikeTiming];
temp = [temp.PeriDP];
S_DP5 = [temp.SpectralPower];

figure('Position', [0 0 300 300])
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
      uc=[0.5 0.5 0 .5];
  end
  shadedErrorBar(log10(f_DP5), ...
    mean(S_DP5(:,u2plot),2), std(S_DP5(:,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', uc, 'LineWidth', 2.5}, 1)
  text(xl(2), yl(2)-.1*range(yl)*(spi-1), utitle, 'FontWeight', 'bold', 'FontSize', fs, 'Color', uc, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
end
set(gca, 'XTick', xt, 'XTickLabel', ft, 'FontSize', fs) % 'YTick', [],
xlim(xl)
ylim(yl)
xlabel(xlab, 'FontSize', fs)
ylabel(ylab, 'FontSize', fs)
title('          DP -150~100ms', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig7Di'))

%% S_DP5 stick IQR VERTICAL
pltopt='loglin';
figure('Position', [200 200 600 300])
for whichband=1:4
  
  switch whichband
    case 1
      fband=f_DP5>=8 & f_DP5<15;
      xlab='8-15 Hz';
      bandtit='Alpha';
    case 2
      fband=f_DP5>=15 & f_DP5<30;
      xlab='15-30 Hz';
      bandtit='Beta';
    case 3
      fband=f_DP5>=30 & f_DP5<=55;
      xlab='30-55 Hz';
      bandtit='Gamma';
    case 4
      fband=f_DP5>=66 & f_DP5<120;
      xlab='66-120 Hz';
      bandtit='High-Gamma';
  end
  
  switch pltopt
    case 'linlin'
      tempS=geomean(S_DP5(fband,:), 1);
      yl=[3*10^-7 3.5*10^-6];
    case 'loglin'
      tempS=geomean(S_DP5(fband,:), 1);
      yl=[3*10^-7 3.5*10^-6];
    case 'loglog'
      tempS=log10(geomean(S_DP5(fband,:), 1));
      yl=log10([5*10^-7 4*10^-6]);
  end
  
  if whichband==1
    yl=[.9*10^-6 7*10^-6];
  end
  
  tempS=10^6*tempS;
  yl=10^6*yl;
  
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
    tempcorr = [tempcorr tempS(u2plot)];
    templab = [templab spi*ones(1,numel(u2plot))];
    
    plot(spi, median(tempS(u2plot)), 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 14)
    plot([spi spi], [prctile(tempS(u2plot),25) prctile(tempS(u2plot),75)], '-', 'Color', uc, 'LineWidth', 4)
    
  end
  [p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
  disp(xlab); disp(p)
  c=multcompare(stats,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), yl(2)-range(yl)/34*[ci+1 ci+1], 'k-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), yl(2)-range(yl)/34*[ci ci], 'k*', 'LineWidth', 1)
  end
  set(gca, 'XTick', 1:spi, 'XTickLabel', xlabs, 'XTickLabelRotation', 45, 'FontSize', fs)
  xlim([0.5 spi+.5])
  ylim(yl)
  ylabel(xlab, 'FontSize', fs)
  title(bandtit, 'FontSize', fs)
  
  text(0.42, yl(2)-.12*range(yl), 'X10^-^6', 'FontAngle', 'italic', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
end
print('-painters', '-dpdf', strcat(figPath, 'fig7Dii'))
