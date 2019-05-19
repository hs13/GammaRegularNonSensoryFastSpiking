%% CV PSTH
temp = [FS.ISI];
temp = [temp.Sliding250ms];
temp = [temp.MatchedHit];
hCVmat = [temp.CV];

temp = [FS.ISI];
temp = [temp.Sliding250ms];
temp = [temp.MatchedMiss];
mCVmat = [temp.CV];

lw=4;
yl=[0.7 1];
figure('Position', [100 100 900 300])
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
  plot([-125 -125], [0 1], '--', 'Color', [0 0 0 0.7])
  shadedErrorBar(FS(1).ISI.param.Sliding250msTimeline, nanmean(mCVmat(:,u2plot),2), ...
    nanstd(mCVmat(:,u2plot),0,2)/sqrt(numel(u2plot)), {'r', 'LineWidth', lw},1)
  shadedErrorBar(FS(1).ISI.param.Sliding250msTimeline, nanmean(hCVmat(:,u2plot),2), ...
    nanstd(hCVmat(:,u2plot),0,2)/sqrt(numel(u2plot)), {'b', 'LineWidth', lw},1)
  set(gca, 'XTick', [-250 -125 0 100 500], 'FontSize', fs)
  xlim([-250 100])
  if spi==1
    ylabel('CV', 'FontSize', fs)
    text(100,yl(1)+.125*range(yl), 'Hit', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'blue', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    text(100,yl(1), 'Miss', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
  end
  xlabel(' ', 'FontSize', fs)
  % if spi==2
  % xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
  % end
  ylim(yl)
  title(strcat(utitle,' (N=',num2str(numel(u2plot)), ')'), 'FontSize', fs)
  if spi==3
    text(110, mean(yl)+.03*range(yl), 'ISI Regularity', 'Rotation', 270, 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    text(110, yl(2)+.05*range(yl), 'low', 'FontAngle', 'italic', 'Rotation', 270, 'FontSize', fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    text(110, yl(1)-.03*range(yl), 'high', 'FontAngle', 'italic', 'Rotation', 270, 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
  end
end
print('-painters', '-dpdf', strcat(figPath, 'fig6A'))

%% logscale inverse ISI colorplot trial normalized probability
temp = [FS.ISI];
temp = [temp.Sliding250ms];
temp = [temp.MatchedHit];
hISImat = cat(3,temp.ISIrate);

temp = [FS.ISI];
temp = [temp.Sliding250ms];
temp = [temp.MatchedMiss];
mISImat = cat(3,temp.ISIrate);

xl=[-250 100];
ft=[8 15 30 60 120];

isihm=hISImat-mISImat;
figure('Position', [100 100 900 300])
% figure('Position', [0 0 800 800])
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
  h=pcolor(FS(1).ISI.param.Sliding250msTimeline, log10(1000./FS(1).ISI.param.NormProbtimeline), ...
    squeeze(nanmean(isihm(:,:,u2plot),3)));
  set(h, 'EdgeColor', 'none');
  hold on
  %   plot([100 100], [0 100], 'w-', 'LineWidth', 1.5)
  xlim(xl)
  ylim([log10(ft(1)) log10(ft(end))])
  %   caxis([-0.075 0.075])
  caxis([-0.04 0.04])
  % %     colorbar
  set(gca, 'XTick', [-250 -125 0 100 500], 'YTick', log10(ft), 'YTickLabel', ft, 'FontSize', fs)
  title(' ', 'FontSize', fs)
  xlabel(' ', 'FontSize', fs)
  title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  if spi==2
    xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
  end
  if spi==1
    ylabel('log ISI^-^1 (s^-^1)', 'FontSize', fs)
  end
end
print('-painters', '-dpdf', strcat(figPath, 'fig6B'))

%% spiking fft spectrogram
temp = [FS.SpikeSpectralPower];
temp = [temp.Sliding250ms];
hSmat = cat(3,temp.MatchedHit);

temp = [FS.SpikeSpectralPower];
temp = [temp.Sliding250ms];
mSmat = cat(3,temp.MatchedMiss);

xl=[-250 100];
ft=[8 15 30 60 120];
figure('Position', [0 0 900 300])
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
  h=pcolor(FS(1).SpikeSpectralPower.param.Sliding250msTimeline, log10(FS(1).SpikeSpectralPower.param.frequency.Sliding250ms), ...
    squeeze(nanmean(hSmat(:,:,u2plot)-mSmat(:,:,u2plot), 3))');
  set(h, 'EdgeColor', 'none');
  xlim(xl)
  ylim([log10(ft(1)) log10(ft(end))])
  set(gca, 'XTick', [-250:125:0 100], 'YTick', log10(ft), 'YTickLabel', ft, 'XGrid', 'on', 'FontSize', fs)
  title(' ', 'FontSize', fs)
  xlabel(' ', 'FontSize', fs)
  %   title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  if spi==2
    xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
  end
  if spi==1
    ylabel('log Frequency (Hz)', 'FontSize', fs)
  end
%   colorbar
  caxis([-1.1 1.1])
end
% colormap jet
print('-painters', '-dpdf', strcat(figPath, 'fig6C'))
