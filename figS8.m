hcol=[1 .5 0];
hlab='False Alarm';
hfield='FalseAlarm';
mcol=[.3 .3 .3];
mlab='Correct Rejection';
mfield='CorrectRejection';

%% PSTH
gausssigma=3;
gausswindow=21;
tempfil= exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter=tempfil/sum(tempfil);

temp = [FS.PSTH];
hpsth = [temp.(hfield)];
mpsth = [temp.(mfield)];

fs=16;
figure('Position', [0 100 1300 300])
for spi=1:3
  yl=[4 25];
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
    case 2
      u2plot=grnsFS;
      utitle='grnsFS';
    case 3
      u2plot=plnsFS;
      utitle='plnsFS';
  end
  
  subplot(1,12,4*(spi-1)+[1:3])
  hold all
  plot(FS(1).PSTH.timeline, nanmean(conv2(hpsth(:,u2plot), tempfilter, 'same'),2),'Color', hcol, 'LineWidth', 2.5)
  plot(FS(1).PSTH.timeline, nanmean(conv2(mpsth(:,u2plot), tempfilter, 'same'),2),'Color', mcol, 'LineWidth', 2.5)
  if spi==1
    text(-250,yl(1)+.1*range(yl), strcat({' '}, hlab), 'FontSize', fs, 'FontWeight', 'bold', 'Color', hcol, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    text(-250,yl(1), strcat({' '}, mlab), 'FontSize', fs, 'FontWeight', 'bold', 'Color', mcol, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
  end
  xlim([-250 100])
  ylim(yl)
  set(gca, 'XTick', [-250 -125 0 100 500], 'FontSize', fs, 'ActivePositionProperty', 'OuterPosition')
  xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
  if spi==1
    ylabel('Rate (Hz)', 'FontSize', fs)
  else
    ylabel(' ', 'FontSize', fs)
  end
  title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  
  
  subplot(1,12,4*spi)
  hold all
  xt=[];
  for tpi=1:2
    switch tpi
      case 1
        hpos=1;
        mpos=2;
        ssctoi=TSTtimeline>=-250 & TSTtimeline<=0;
      case 2
        hpos=3.5;
        mpos=4.5;
        ssctoi=TSTtimeline>=0 & TSTtimeline<=100;
    end
    hp=mean(hpsth(ssctoi,:),1);
    mp=mean(mpsth(ssctoi,:),1);
    errorbar(hpos, nanmedian(hp(u2plot)), prctile(hp(u2plot),50)-prctile(hp(u2plot),25), prctile(hp(u2plot),75)-prctile(hp(u2plot),50), ...
      's-', 'Color', hcol, 'LineWidth', 2.5, 'MarkerFaceColor', hcol, 'MarkerSize',10, 'CapSize', 0)
    errorbar(mpos, nanmedian(mp(u2plot)), prctile(mp(u2plot),50)-prctile(mp(u2plot),25), prctile(mp(u2plot),75)-prctile(mp(u2plot),50), ...
      's-', 'Color', mcol, 'LineWidth', 2.5, 'MarkerFaceColor', mcol, 'MarkerSize',10, 'CapSize', 0)
    
    pval=signrank(hp(u2plot), mp(u2plot));
    if pval<0.05/3
      plot((hpos+mpos)/2, yl(2), 'k*', 'MarkerSize', 5,'LineWidth',1)
    end
    xt=[xt (hpos+mpos)/2];
  end
  xlim([.5 5])
  ylim(yl)
  set(gca, 'YTick', [], 'XTick', [1.5 4], 'XTickLabel', {' ', ' '}, 'XTickLabelRotation', 45, 'FontSize', fs, 'ActivePositionProperty', 'OuterPosition')
  hpre=text(1.5, yl(1), 'Pre', 'FontSize', 14, 'VerticalAlignment', 'top');
  set(hpre, 'Rotation', -45)
  hpost=text(4, yl(1), 'Post', 'FontSize', 14, 'VerticalAlignment', 'top');
  set(hpost, 'Rotation', -45)
  
  xlabel(' ', 'FontSize', fs)
  title(' ', 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS8psth'))

%% LFP PSD
fs=16;
for pltopt=1:2
switch pltopt
  case 1
    ft=[8 15 30 60 120];
    yl=[-.2 2.5];
    xf=log10(Session(1).LFPPower.param.frequency.Pre250ms);
    temp = [Session.LFPPower];
    temp = [temp.Pre250ms];
    xtit = 'Aileft';
  case 2
    ft=[20 30 60 120];
    yl=[0 2.5];
    xf=log10(Session(1).LFPPower.param.frequency.Post100ms);
    temp = [Session.LFPPower];
    temp = [temp.Post100ms];
    xtit = 'Bileft';
end

savgh=log10([temp.(hfield)]);
savgm=log10([temp.(mfield)]);

figure('Position',[0 0 300 300])
hold all
shadedErrorBar(xf, nanmean(savgh,2), nanstd(savgh,0,2)/sqrt(size(savgh,2)), {'Color', hcol, 'LineWidth', 2.5}, 1)
shadedErrorBar(xf, nanmean(savgm,2), nanstd(savgm,0,2)/sqrt(size(savgm,2)), {'Color', mcol, 'LineWidth', 2.5}, 1)
set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs)
xlim(log10([ft(1) ft(end)]))
ylim(yl)
title(strcat('LFP (N=', num2str(size(savgh,2)), ')'))
xlabel('log Frequency (Hz)', 'FontSize', fs)
ylabel('log Power', 'FontSize', fs)
% if pltopt==1
%   text(log10(ft(end)),yl(2), hlab, 'Color', hcol, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
%   text(log10(ft(end)),yl(2)-.1*range(yl), mlab, 'Color', mcol, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
% end
print('-painters', '-dpdf', strcat(figPath, 'figS8', xtit))
end

%% LFP alpha beta gamma hits & misses
for pltopt=1:2
switch pltopt
  case 1
    ft=[8 15 30 60 120];
    yl=[-.2 2.5];
    xf=Session(1).LFPPower.param.frequency.Pre250ms;
    temp = [Session.LFPPower];
    temp = [temp.Pre250ms];
    xtit = 'Airight';
    ylab='Pre 250ms';
  case 2
    ft=[20 30 60 120];
    yl=[0 2.5];
    xf=Session(1).LFPPower.param.frequency.Post100ms;
    temp = [Session.LFPPower];
    temp = [temp.Post100ms];
    xtit = 'Biright';
    ylab='Post 100ms';
end

savgh=[temp.(hfield)];
savgm=[temp.(mfield)];

figure('Position',[100 100 200 300])
ut={};
hold all
for spi=1:3
  switch spi
    case 3
      foi=xf>=66 & xf<=120;
      rtitle='66~120Hz';
    case 2
      foi=xf>=30 & xf<=55;
      rtitle='30~55Hz';
    case 1
      if pltopt==2
        foi=xf>=20 & xf<30;
        rtitle='20~30Hz';
      else
        foi=xf>=15 & xf<30;
        rtitle='15~30Hz';
      end
  end
  hp=log10(geomean(savgh(foi,:), 1));
  mp=log10(geomean(savgm(foi,:), 1));
  
  errorbar(2.5*spi-1, nanmedian(hp), prctile(hp,50)-prctile(hp,25), prctile(hp,75)-prctile(hp,50), ...
    's-', 'Color', hcol, 'MarkerFaceColor', hcol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  errorbar(2.5*spi, nanmedian(mp), prctile(mp,50)-prctile(mp,25), prctile(mp,75)-prctile(mp,50), ...
    's-', 'Color', mcol, 'MarkerFaceColor', mcol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  pval=signrank(hp, mp);
  if pval<0.05
    plot(2.5*spi-0.5,yl(2), 'k*', 'LineWidth', 1)
  end
  ut=[ut rtitle];
end
ylim(yl)
xlim([1 2.5*spi+0.5])
set(gca, 'XTick', 2.5*[1:3]-0.5, 'XTickLabel', ut, 'XTickLabelRotation', 45, 'FontSize', fs)
ylabel('log Power', 'FontSize', fs)
%   xlabel(' ', 'FontSize', fs)
title(ylab, 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS8', xtit))
end

%% stick IQR CV vs meanISI vs stdISI hits & misses
for pltopt=1:2
switch pltopt
  case 1
    ylab='Pre 250ms';
    temp = [FS.ISI];
    temp = [temp.Pre250ms];
    temp = [temp.(hfield)];
    hp=[temp.CV];
    temp = [FS.ISI];
    temp = [temp.Pre250ms];
    temp = [temp.(mfield)];
    mp=[temp.CV];
    xtit = 'Aii';
    yl=[.5 1];
  case 2
    ylab='Post 100ms';
    temp = [FS.ISI];
    temp = [temp.Post100ms_2ndST];
    temp = [temp.(hfield)];
    hp=[temp.CV];
    temp = [FS.ISI];
    temp = [temp.Post100ms_2ndST];
    temp = [temp.(mfield)];
    mp=[temp.CV];
    xtit = 'Bii';
    yl=[0 1.4];
end
rtitle='CV';

figure('Position',[100 100 200 300])
ut={};
hold all
for spi=1:3
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
    case 2
      u2plot=grnsFS;
      utitle='grnsFS';
    case 3
      u2plot=plnsFS;
      utitle='plnsFS';
  end
  errorbar(2.5*spi-1, nanmedian(hp(u2plot)), prctile(hp(u2plot),50)-prctile(hp(u2plot),25), prctile(hp(u2plot),75)-prctile(hp(u2plot),50), ...
    's-', 'Color', hcol, 'MarkerFaceColor', hcol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  errorbar(2.5*spi, nanmedian(mp(u2plot)), prctile(mp(u2plot),50)-prctile(mp(u2plot),25), prctile(mp(u2plot),75)-prctile(mp(u2plot),50), ...
    's-', 'Color', mcol, 'MarkerFaceColor', mcol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  %    effect size
  Scohen=sqrt((numel(u2plot)*var(hp(u2plot))+numel(u2plot)*var(mp(u2plot)))/(2*numel(u2plot)));
  Dcohen=(mean(hp(u2plot))-mean(mp(u2plot)))/Scohen;
  pval=signrank(hp(u2plot), mp(u2plot));
  if pval<0.05/3
    plot(2.5*spi-0.5,yl(2), 'k*', 'LineWidth', 1)
  end
  ut=[ut utitle];
end
ylim(yl)
xlim([1 2.5*spi+0.5])
set(gca, 'XTick', 2.5*[1:3]-0.5, 'XTickLabel', ut, 'XTickLabelRotation', 45, 'FontSize', fs)
ylabel(rtitle, 'FontSize', fs)
title(ylab, 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS8', xtit))
end

%% prestim ISI histogram hits & misses figure
for pltopt=1:2
  switch pltopt
    case 1
      ft=[8 15 30 60 120];
      temp = [FS.ISI];
      temp = [temp.Pre250ms];
      temp = [temp.(hfield)];
      hmat=[temp.ISIrate];
      temp = [FS.ISI];
      temp = [temp.Pre250ms];
      temp = [temp.(mfield)];
      mmat=[temp.ISIrate];
      xtit = 'Aiii';
    case 2
      ft=[20 30 60 120];
      temp = [FS.ISI];
      temp = [temp.Post100ms_2ndST];
      temp = [temp.(hfield)];
      hmat=[temp.ISIrate];
      temp = [FS.ISI];
      temp = [temp.Post100ms_2ndST];
      temp = [temp.(mfield)];
      mmat=[temp.ISIrate];
      xtit = 'Biii';
  end
  
  xisi=log10(1000./FS(1).ISI.param.NormProbtimeline);
  xl=log10([ft(1) ft(end)]);
  
  figure('Position',[100 100 900 300])
  for spi=1:3
    switch spi
      case 1
        u2plot=sFS;
        utitle='sFS';
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
    end
    if pltopt==1
      yl=[0 0.27];
    yt=[0 .25];
    else
      yl=[0 0.5];
    yt=yl;
    end
    
    subplot(1,3, spi)
    hold all
    shadedErrorBar(xisi, nanmean(mmat(:,u2plot),2), nanstd(mmat(:,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', mcol, 'LineWidth', 2.5}, 1)
    shadedErrorBar(xisi, nanmean(hmat(:,u2plot),2), nanstd(hmat(:,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', hcol, 'LineWidth', 2.5}, 1)
    xlim(xl)
    ylim(yl)
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', yt, 'FontSize', fs)
    if spi==2
      xlabel('log ISI^-^1 (s^-^1)', 'FontSize', fs)
    end
    if spi==1
      ylabel('Normalized Rate', 'FontSize', fs)
    end
    title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  end
  orient(gcf, 'landscape')
  print('-painters', '-dpdf', strcat(figPath, 'figS8', xtit))
end

%% S2 PSD hit vs miss
for pltopt=1:2
  switch pltopt
    case 1
      ft=[8 15 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Pre250ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Pre250ms];
      hmat = [temp.(hfield)];
      mmat = [temp.(mfield)];
      xtit = 'Avi';
    case 2
      ft=[20 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Post100ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Post100ms];
      hmat = [temp.(hfield)];
      mmat = [temp.(mfield)];
      xtit = 'Bvi';
  end
      yl=[6 16];
  
  
  figure('Position',[0 0 900 300])
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
    subplot(1,3,spi)
    hold all
    lw=2.5;
    shadedErrorBar(log10(fVec), (nanmean(mmat(:,u2plot),2)), (nanstd(mmat(:,u2plot),0,2))/sqrt(numel(u2plot)), {'Color', mcol, 'LineWidth', lw}, 1)
    shadedErrorBar(log10(fVec), (nanmean(hmat(:,u2plot),2)), (nanstd(hmat(:,u2plot),0,2))/sqrt(numel(u2plot)), {'Color', hcol, 'LineWidth', lw}, 1)
    ylim(yl)
    xlim([log10(ft(1)) log10(ft(end))])
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', 0:5:20, 'FontSize', fs) % , 'YTick', yl
    if spi==1
      % text(log10(ft(end)),yl(2), hlab, 'Color', hcol, 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
      % text(log10(ft(end)),yl(2)-.1*range(yl), mlab, 'Color', mcol, 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
      ylabel('Spike Spectral Power', 'FontSize', fs)
    end
    if spi==2
      xlabel('log Frequency (Hz)', 'FontSize', fs)
    end
    title(' ', 'FontSize', fs)
    % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  end
  orient(gcf, 'landscape')
  print('-painters', '-dpdf', strcat(figPath, 'figS8', xtit))
end
