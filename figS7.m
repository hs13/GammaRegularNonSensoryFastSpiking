hcol=[0 0 1];
hlab='Hit';
hfield='MatchedHit';
mcol=[1 0 0];
mlab='Miss';
mfield='MatchedMiss';

%% LFP PSD
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
if pltopt==1
  text(log10(ft(end)),yl(2), hlab, 'Color', hcol, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
  text(log10(ft(end)),yl(2)-.1*range(yl), mlab, 'Color', mcol, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
end
print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
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
  %    effect size
  Scohen=sqrt((length(hp)*var(hp)+length(mp)*var(mp))/(length(hp)+length(mp)));
  Dcohen=(mean(hp)-mean(mp))/Scohen;
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
print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
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
    yl=[.7 1];
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
    yl=[.8 1.6];
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
print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
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
        if pltopt==1
          yl = [0 0.33];
          yt = [0 0.3];
        else
          yl = [0 0.55];
          yt = [0 0.5];
        end
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
        if pltopt==1
          yl = [0 0.2];
          yt=yl;
        else
          yl = [0 0.25];
          yt=yl;
        end
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
        if pltopt==1
          yl = [0 0.2];
          yt=yl;
        else
          yl = [0 0.25];
          yt=yl;
        end
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
  print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
end

%% stick IQR ISI hits & misses
for pltopt=1:2
  
  figure('Position',[100 200*pltopt 900 200])
  
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
    ut={};
    subplot(1,3,spi)
    hold all
    
    for ii=1:3
      if pltopt==1
        temp = [FS.ISI];
        temp = [temp.Pre250ms];
        htemp = [temp.(hfield)];
        mtemp = [temp.(mfield)];
        
        yl=[0 7.2];
        ylab='Pre 250ms ISI Norm Rate';
        xtit = 'Aiv';
      elseif pltopt==2
        temp = [FS.ISI];
        temp = [temp.Post100ms_2ndST];
        htemp = [temp.(hfield)];
        mtemp = [temp.(mfield)];
        
        yl=[0 11];
        ylab='Post 100ms ISI Norm Rate';
        xtit = 'Biv';
      else
        error('specify pltopt')
      end
      hcol=[0 0 1];
      mcol=[1 0 0];
      
      switch ii
        case 1
          hp=[htemp.shorterISIrate];
          mp=[mtemp.shorterISIrate];
          rtitle='<18ms';
        case 2
          hp=[htemp.gammaISIrate];
          mp=[mtemp.gammaISIrate];
          rtitle='18~33ms';
        case 3
          hp=[htemp.longerISIrate];
          mp=[mtemp.longerISIrate];
          rtitle='>33ms';
      end
      
      plot(nanmedian(hp(u2plot)), 2.5*ii-1, 's', 'Color', hcol, 'MarkerFaceColor', hcol, 'MarkerSize', 10)
      plot([prctile(hp(u2plot),25) prctile(hp(u2plot),75)], 2.5*[ii ii]-1, '-', 'Color', hcol, 'LineWidth', 2.5)
      plot(nanmedian(mp(u2plot)), 2.5*ii, 's', 'Color', mcol, 'MarkerFaceColor', mcol, 'MarkerSize', 10)
      plot([prctile(mp(u2plot),25) prctile(mp(u2plot),75)], 2.5*[ii ii], '-', 'Color', mcol, 'LineWidth', 2.5)
      pval=signrank(hp(u2plot), mp(u2plot));
      if pval<0.05/3
        plot(yl(2),2.5*ii-0.5, 'k*', 'LineWidth', 1)
      end
      ut=[ut rtitle];
    end
    if spi==1
      set(gca, 'YTick', 2.5*[1:3]-0.5, 'YTickLabel', ut, 'FontSize', fs)
    else
      set(gca, 'YTick', [], 'FontSize', fs)
    end
    xlim(yl)
    ylim([.75 2.5*ii+0.75])
    if spi==2
      xlabel(ylab, 'FontSize', fs)
    end
    title(' ', 'FontSize', fs)
    %     title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  end
  orient(gcf, 'landscape')
  print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
end

%% shuffle CI % significant units
for pltopt=1:2
  if pltopt==1
    temp = [FS.ISI];
    temp = [temp.Pre250ms];
    xtit = 'Av';
    ft=[8 15 30 60 120];
  elseif pltopt==2
    temp = [FS.ISI];
    temp = [temp.Post100ms_2ndST];
    xtit = 'Bv';
    ft=[20 30 60 120];
  else
    error('specify pltopt')
  end
  htemp = [temp.MatchedHit];
  mtemp = [temp.MatchedMiss];
  hmdiff = [htemp.ISIrate]-[mtemp.ISIrate];
  
  citemp = [temp.HitMissDiff_95CI];
  hmdiffCI = cat(3,citemp.ISIrate);
  sighit = hmdiff > squeeze(hmdiffCI(:,2,:));
  sigmiss = hmdiff < squeeze(hmdiffCI(:,1,:));
  
  xisi=log10(1000./(FS(1).ISI.param.NormProbtimeline-.5));
  xl=log10([ft(1) ft(end)]);
  
  figure('Position', [100 300*(pltopt-1) 900 300])
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
    subplot(1,3,spi)
    hold all
    stairs(xisi, 100*sum(sighit(:,u2plot),2)/numel(u2plot), 'b-', 'LineWidth', 1.5)
    stairs(xisi, 100*sum(sigmiss(:,u2plot),2)/numel(u2plot), 'r-', 'LineWidth', 1.5)
    plot([xisi(1) xisi(end)], [2.5 2.5], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
    if spi<3
      yl=ylim;
    end
    yl(2)=max(yl(2),7);
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs)
    if spi==2
      xlabel('log ISI^-^1 (s^-^1)', 'FontSize', fs)
    end
    if spi==1
      text(xl(1),yl(2), ' Hit>Miss', 'Color', 'blue', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
      text(xl(1),yl(2)-.1*range(yl), ' Miss>Hit', 'Color', 'red', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
      ylabel('% Significant Units', 'FontSize', fs)
    end
    xlim(xl)
    ylim(yl)
    title(' ', 'FontSize', fs)
    % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'))
  end
  orient(gcf, 'landscape')
  print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
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
      yl=[6 16];
    case 2
      ft=[20 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Post100ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Post100ms];
      hmat = [temp.(hfield)];
      mmat = [temp.(mfield)];
      xtit = 'Bvi';
      yl=[7 19];
  end
  
  
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
  print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
end

%% stick IQR spike spectral power hits & misses
for pltopt=1:2
  switch pltopt
    case 1
      ft=[8 15 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Pre250ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Pre250ms];
      hmat = [temp.(hfield)];
      mmat = [temp.(mfield)];
      xtit = 'Avii';
      ylab='Pre 250ms Power';
      yl=[3 20];
    case 2
      ft=[20 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Post100ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Post100ms];
      hmat = [temp.(hfield)];
      mmat = [temp.(mfield)];
      xtit = 'Bvii';
      ylab='Post 100ms Power';
      yl=[3 24];
  end
  
  figure('Position',[100 200*pltopt 900 200])
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
    subplot(1,3,spi)
    ut={};
    hold all
    for ii=1:3
      switch ii
        case 1
          foi=fVec>=66 & fVec<=120;
          rtitle='66~120Hz';
        case 2
          foi=fVec>=30 & fVec<=55;
          rtitle='30~55Hz';
        case 3
          if contains(ylab, '100ms')
            foi=fVec>=20 & fVec<30;
            rtitle='20~30Hz';
          else
            foi=fVec>=15 & fVec<30;
            rtitle='15~30Hz';
          end
      end
      
      hp=squeeze(geomean(hmat(foi,:),1));
      hcol=[0 0 1];
      mp=squeeze(geomean(mmat(foi,:),1));
      mcol=[1 0 0];
      
      plot(nanmedian(hp(u2plot)), 2.5*ii-1, 's', 'Color', hcol, 'MarkerFaceColor', hcol, 'MarkerSize', 10)
      plot([prctile(hp(u2plot),25) prctile(hp(u2plot),75)], 2.5*[ii ii]-1, '-', 'Color', hcol, 'LineWidth', 2.5)
      plot(nanmedian(mp(u2plot)), 2.5*ii, 's', 'Color', mcol, 'MarkerFaceColor', mcol, 'MarkerSize', 10)
      plot([prctile(mp(u2plot),25) prctile(mp(u2plot),75)], 2.5*[ii ii], '-', 'Color', mcol, 'LineWidth', 2.5)
      pval=signrank(hp(u2plot), mp(u2plot));
      if pval<0.05/3
        plot(yl(2),2.5*ii-0.5, 'k*', 'LineWidth', 1)
      end
      ut=[ut rtitle];
    end
    if spi==1
      set(gca, 'YTick', 2.5*[1:3]-0.5, 'YTickLabel', ut, 'FontSize', fs)
    else
      set(gca, 'YTick', [], 'FontSize', fs)
    end
    xlim(yl)
    ylim([.75 2.5*ii+0.75])
    if spi==2
      xlabel(ylab, 'FontSize', fs)
    end
    title(' ', 'FontSize', fs)
    %   title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  end
  orient(gcf, 'landscape')
  print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
end

%% shuffle CI % significant units
for pltopt=1:2
  switch pltopt
    case 1
      ft=[8 15 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Pre250ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Pre250ms];
      xtit = 'Aviii';
    case 2
      ft=[20 30 60 120];
      fVec = FS(1).SpikeSpectralPower.param.frequency.Post100ms;
      temp = [FS.SpikeSpectralPower];
      temp = [temp.Post100ms];
      xtit = 'Bviii';
  end
  hmdiff = [temp.MatchedHit] - [temp.MatchedMiss];
  hmdiffCI = cat(3,temp.HitMissDiff_95CI);
  sighit = hmdiff > squeeze(hmdiffCI(:,2,:));
  sigmiss = hmdiff < squeeze(hmdiffCI(:,1,:));
  
  figure('Position', [100 300*(pltopt-1) 900 300])
  for spi=1:3
    switch spi
      case 1
        u2plot=sFS;
        utitle='sFS';
        yl=[0 42];
      case 2
        u2plot=grnsFS;
        utitle='grnsFS';
        yl=[0 20];
      case 3
        u2plot=plnsFS;
        utitle='plnsFS';
        yl=[0 20];
    end
    subplot(1,3,spi)
    hold all
    stairs(log10(fVec), 100*sum(sigmiss(:,u2plot),2)/numel(u2plot), '-', 'Color', [1 0 0 .5], 'LineWidth', 1.5)
    stairs(log10(fVec), 100*sum(sighit(:,u2plot),2)/numel(u2plot), '-', 'Color', [0 0 1 .5], 'LineWidth', 1.5)
    plot([0 250], [2.5 2.5], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
    % yl=ylim;
    if spi==1
      % text(log10(ft(end)),yl(2), 'Hit>Miss', 'Color', 'blue', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
      % text(log10(ft(end)),yl(2)-.1*range(yl), 'Miss>Hit', 'Color', 'red', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
      ylabel('% Significant Units', 'FontSize', fs)
    end
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs)
    ylim(yl)
    xlabel('Frequency (Hz)', 'FontSize', fs)
    title(' ', 'FontSize', fs)
    % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'))
    xlim([log10(ft(1)) log10(ft(end))])
  end
  orient(gcf, 'landscape')
  print('-painters', '-dpdf', strcat(figPath, 'figS7', xtit))
end