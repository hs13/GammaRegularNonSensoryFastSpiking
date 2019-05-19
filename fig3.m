load(strcat(dataPath, 'Session.mat'))

%% LFP PSD
yl=[-.2 2.5];

ft=[20 30 60 120];
fvec=Session(1).LFPPower.param.frequency.Post100ms;
xf=log10(fvec);

temp = [Session.LFPPower];
temp = [temp.Post100ms];
savgpost = [temp.Maxstim];

temp = [Session.LFPPower];
temp = [temp.Pre100ms];
savgpre = [temp.Maxstim];

postcol=[1 0 1];
precol=[0 0 0];
postlab='Post';
prelab='Pre';

figure('Position',[100 100 300 300])
hold all
shadedErrorBar(xf, nanmean(log10(savgpost),2), nanstd(log10(savgpost),0,2)/sqrt(size(savgpost,2)), {'Color', postcol, 'LineWidth', 2.5}, 1)
shadedErrorBar(xf, nanmean(log10(savgpre),2), nanstd(log10(savgpre),0,2)/sqrt(size(savgpre,2)), {'Color', precol, 'LineWidth', 2.5}, 1)
set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', 0:2, 'FontSize', fs)
xlim(log10([ft(1) ft(end)]))
ylim(yl)
text(log10(ft(end)),yl(2)-.1*range(yl), postlab, 'Color', postcol, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
text(log10(ft(end)),yl(2), prelab, 'Color', precol, 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
title(strcat('LFP (N=', num2str(size(savgpost,2)), ')'))
xlabel('log Frequency (Hz)', 'FontSize', fs)
ylabel('log Power', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig3Ai'))


%% stick IQR LFP power hits & misses
yl=[-.2 2.5];

figure('Position',[100 100 200 300])
ut={};
hold all
for spi=1:3
  switch spi
    case 3
      foi=fvec>=66 & fvec<=120;
      rtitle='66~120Hz';
    case 2
      foi=fvec>=30 & fvec<=55;
      rtitle='30~55Hz';
    case 1
      if isequal(fvec, Session(1).LFPPower.param.frequency.Post100ms)
        foi=fvec>=20 & fvec<30;
        rtitle='20~30Hz';
      else
        foi=fvec>=15 & fvec<30;
        rtitle='15~30Hz';
      end
  end
  postvec=log10(geomean(savgpost(foi,:), 1));
  prevec=log10(geomean(savgpre(foi,:), 1));
  
  errorbar(2.5*spi-1, nanmedian(postvec), prctile(postvec,50)-prctile(postvec,25), prctile(postvec,75)-prctile(postvec,50), ...
    's-', 'Color', postcol, 'MarkerFaceColor', postcol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  errorbar(2.5*spi, nanmedian(prevec), prctile(prevec,50)-prctile(prevec,25), prctile(prevec,75)-prctile(prevec,50), ...
    's-', 'Color', precol, 'MarkerFaceColor', precol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  %    effect size
  Scohen=sqrt((length(postvec)*var(postvec)+length(prevec)*var(prevec))/(length(postvec)+length(prevec)));
  Dcohen=(mean(postvec)-mean(prevec))/Scohen;
  pval=signrank(postvec, prevec);
  if pval<0.05/3
    plot(2.5*spi-0.5,yl(2), 'k*', 'LineWidth', 1)
  end
  ut=[ut rtitle];
end
ylim(yl)
xlim([1 2.5*spi+0.5])
set(gca, 'XTick', 2.5*[1:3]-0.5, 'XTickLabel', ut, 'XTickLabelRotation', 45, 'FontSize', fs)
ylabel('log Power', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig3Aii'))

%% stick IQR CV vs meanISI vs stdISI pre & post
yl=[.8 1.7];
rtitle='CV';
temp = [FS.ISI];
temp = [temp.Post100ms_2ndST];
temp = [temp.Maxstim];
postvec = [temp.CV];

temp = [FS.ISI];
temp = [temp.Pre100ms_2ndST];
temp = [temp.Maxstim];
prevec = [temp.CV];

postcol=[1 0 1];
precol=[0 0 0];

%   subplot(1,3,ii)
figure('Position',[100 100 200 300])
ut={};
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
  errorbar(2.5*spi-1, nanmedian(postvec(u2plot)), prctile(postvec(u2plot),50)-prctile(postvec(u2plot),25), prctile(postvec(u2plot),75)-prctile(postvec(u2plot),50), ...
    's-', 'Color', postcol, 'MarkerFaceColor', postcol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  errorbar(2.5*spi, nanmedian(prevec(u2plot)), prctile(prevec(u2plot),50)-prctile(prevec(u2plot),25), prctile(prevec(u2plot),75)-prctile(prevec(u2plot),50), ...
    's-', 'Color', precol, 'MarkerFaceColor', precol, 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2.5)
  %    effect size
  Scohen=sqrt((numel(u2plot)*var(postvec(u2plot))+numel(u2plot)*var(prevec(u2plot)))/(2*numel(u2plot)));
  Dcohen=(mean(postvec(u2plot))-mean(prevec(u2plot)))/Scohen;
  pval=signrank(postvec(u2plot), prevec(u2plot));
  if pval<0.05/3
    plot(2.5*spi-0.5,yl(2), 'k*', 'LineWidth', 1)
  end
  ut=[ut utitle];
end
%     yl=ylim;
ylim(yl)
xlim([1 2.5*spi+0.5])
set(gca, 'XTick', 2.5*[1:3]-0.5, 'XTickLabel', ut, 'XTickLabelRotation', 45, 'FontSize', fs)
ylabel(rtitle, 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig3B'))

%% pre post ISI average across cell group
xisi=log10(1000./FS(1).ISI.param.NormProbtimeline);
ft=[8 15 30 60 120];
xl=log10([ft(1) ft(end)]);

temp = [FS.ISI];
temp = [temp.Post100ms_2ndST];
temp = [temp.Maxstim];
postisi = [temp.ISIrate];

temp = [FS.ISI];
temp = [temp.Pre100ms_2ndST];
temp = [temp.Maxstim];
preisi = [temp.ISIrate];

figure('Position', [0 0 900 300])
for spi=1:3
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
      yl=[0 0.7];
      yt=0:0.2:0.8;
    case 2
      u2plot=grnsFS;
      utitle='grnsFS';
      yl=[0 0.22];
      yt=[0 0.1 0.2];
    case 3
      u2plot=plnsFS;
      utitle='plnsFS';
      yl=[0 0.22];
      yt=[0 0.1 0.2];
  end
  subplot(1,3, spi)
  hold all
  shadedErrorBar(xisi, nanmean(preisi(:,u2plot),2), nanstd(preisi(:,u2plot),0,2)/sqrt(numel(u2plot)), {'k-', 'LineWidth', 2.5}, 1)
  shadedErrorBar(xisi, nanmean(postisi(:,u2plot),2), nanstd(postisi(:,u2plot),0,2)/sqrt(numel(u2plot)), {'m-', 'LineWidth', 2.5}, 1)
  if spi==1
    text(xl(2), yl(1)+.1*range(yl), 'Pre', 'Color', 'k', 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
    text(xl(2), yl(1), ' Post', 'Color', 'm', 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
  end
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
print('-painters', '-dpdf', strcat(figPath, 'fig3C'))

%% pre post ISI shuffle difference % significant units
xisi=log10(1000./FS(1).ISI.param.NormProbtimeline);
ft=[8 15 30 60 120];
xl=log10([ft(1) ft(end)]);

postprediff = postisi-preisi;

temp = [FS.ISI];
temp = [temp.PostPreDiff100ms_2ndST_95CI];
temp = [temp.Maxstim];
postprediffCI = cat(3,temp.ISIrate);
sigpre = postprediff < squeeze(postprediffCI(:,1,:));
sigpost = postprediff > squeeze(postprediffCI(:,2,:));

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
  yl=[0 52];
  subplot(1,3,spi)
  hold all
  stairs(xisi, 100*sum(sigpre(:,u2plot),2)/numel(u2plot), 'k-', 'LineWidth', 1.5)
  stairs(xisi, 100*sum(sigpost(:,u2plot),2)/numel(u2plot), 'm-', 'LineWidth', 1.5)
  plot([xisi(1) xisi(end)], [2.5 2.5], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
  set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs)
  if spi==2
    xlabel('log ISI^-^1 (s^-^1)', 'FontSize', fs)
  end
  if spi==1
    text(xl(1),yl(2), ' Pre>Post', 'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
    text(xl(1),yl(2)-.1*range(yl), ' Post>Pre', 'Color', 'm', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
    ylabel('% Significant Units', 'FontSize', fs)
  end
  xlim(xl)
  ylim(yl)
  title(' ','FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'fig3D'))

%% S2 PSD post vs pre
ft=[20 30 60 120];
fvec=FS(1).SpikeSpectralPower.param.frequency.Post100ms;
xf=log10(fvec);

temp = [FS.SpikeSpectralPower];
temp = [temp.Post100ms];
savgpost = [temp.Maxstim];
postcol=[1 0 1];
postlab='Post';

temp = [FS.SpikeSpectralPower];
temp = [temp.Pre100ms];
savgpre = [temp.Maxstim];
precol=[0 0 0];
prelab='Pre';
yl=[6 22];

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
  shadedErrorBar(xf, (nanmean(savgpre(:,u2plot),2)), (nanstd(savgpre(:,u2plot),0,2))/sqrt(numel(u2plot)), {'Color', precol, 'LineWidth', lw}, 1)
  shadedErrorBar(xf, (nanmean(savgpost(:,u2plot),2)), (nanstd(savgpost(:,u2plot),0,2))/sqrt(numel(u2plot)), {'Color', postcol, 'LineWidth', lw}, 1)
  ylim(yl)
  xlim([log10(ft(1)) log10(ft(end))])
  set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', 0:5:20, 'FontSize', fs) % , 'YTick', yl
  if spi==1
    text(log10(ft(end)),yl(1), postlab, 'Color', postcol, 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    text(log10(ft(end)),yl(1)+.1*range(yl), prelab, 'Color', precol, 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    ylabel('Spike Spectral Power', 'FontSize', fs)
  end
  if spi==2
    xlabel('log Frequency (Hz)', 'FontSize', fs)
  end
  title(' ', 'FontSize', fs)
  % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'fig3E'))

%% pre post ISI shuffle difference % significant units
ft=[20 30 60 120];
fvec=FS(1).SpikeSpectralPower.param.frequency.Post100ms;
xf=log10(fvec);
xl=log10([ft(1) ft(end)]);

postprediff = savgpost-savgpre;

temp = [FS.SpikeSpectralPower];
temp = [temp.PostPreDiff100ms_95CI];
postprediffCI = cat(3,temp.Maxstim);
sigpreS2 = postprediff < squeeze(postprediffCI(:,1,:));
sigpostS2 = postprediff > squeeze(postprediffCI(:,2,:));

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
  yl=[0 100];
  subplot(1,3,spi)
  hold all
  stairs(xf, 100*sum(sigpreS2(:,u2plot),2)/numel(u2plot), 'k-', 'LineWidth', 1.5)
  stairs(xf, 100*sum(sigpostS2(:,u2plot),2)/numel(u2plot), 'm-', 'LineWidth', 1.5)
  plot([xf(2) xf(end)], [2.5 2.5], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
  set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'FontSize', fs)
  if spi==2
    xlabel('log Frequency (Hz)', 'FontSize', fs)
  end
  if spi==1
    text(xl(1),yl(1)+.1*range(yl), ' Pre>Post', 'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    text(xl(1),yl(1), ' Post>Pre', 'Color', 'm', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    ylabel('% Significant Units', 'FontSize', fs)
  end
  xlim(xl)
  ylim(yl)
  title(' ', 'FontSize', fs)
  % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'),'FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'fig3F'))
