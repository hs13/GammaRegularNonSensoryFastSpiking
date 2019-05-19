temp = [FS.FanoFactor];
temp=[temp.Pre1s];
ffvec=[temp.AllTrials];

%% FF stick IQR horizontal + histogram
tempvar = ffvec;
he=0:0.4:10;
xl=[.5 8];
xlab='Fano Factor';

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [0 200 250 300])
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
  plot(xl(2)-range(xl)*0.05*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', .5)
  plot(xl(2)+range(xl)*0.02-range(xl)*0.05*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'LineWidth', .5)
end
axis off
set(gca, 'XTick', [], 'YTick', [], 'FontSize', fs)%, 'ActivePositionProperty', 'outerposition')
ylim([0.5 spi+.5])
xlim(xl)

yl=[0 20];
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
  text(xl(2)-.1*range(xl), yl(2)-.1*range(yl)*(spi-1), utitle, ...
    'Color', uc, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
end
set(gca, 'FontSize', fs, 'ActivePositionProperty', 'outerposition') % 'XTick', .8:.2:1.6,
xlim(xl)
ylim(yl)
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS2A'))

%% scatter FF - CV2
temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
cvsqvec=[temp.CV].^2;

tci=isitcind.pre1s;
ttype='alltrials';

figure('Position', [100 100 300 300])
hold all
plot([0 10], [0 10], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
scatter(ffvec(sFS),cvsqvec(sFS),25, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerFaceAlpha', 0.4)
scatter(ffvec(grnsFS),cvsqvec(grnsFS), 25, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0.4)
scatter(ffvec(plnsFS),cvsqvec(plnsFS), 25, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0], 'MarkerFaceAlpha', 0.4)

text(10,0.5, sprintf('R=%.2f', corr(ffvec',cvsqvec')), ...
  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 16)
ylabel('CV^2', 'FontSize', fs)
xlabel('Fano Factor', 'FontSize', fs)
set(gca, 'FontSize', fs)
axis([0 10 0.5 2.5])
print('-painters', '-dpdf', strcat(figPath, 'figS2Bi'))

%% CVsq/FF stick IQR horizontal + histogram

Rffcvsq=(cvvec.^2)./ffvec;
tempvar=Rffcvsq;
he=0:0.05:2;
xl=[0 1];
% xlab='CV';
xlab='CV^2 / Fano Factor';

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [0 200 250 300])
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
  plot(xl(2)-range(xl)*0.05*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', .5)
  plot(xl(2)+range(xl)*0.02-range(xl)*0.05*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'LineWidth', .5)
end
axis off
set(gca, 'XTick', [], 'YTick', [], 'FontSize', fs)%, 'ActivePositionProperty', 'outerposition')
ylim([0.5 spi+.5])
xlim(xl)

yl=[0 20];
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
  text(xl(2)-.1*range(xl), yl(2)-.1*range(yl)*(spi-1), utitle, ...
    'Color', uc, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
end
set(gca, 'FontSize', fs, 'ActivePositionProperty', 'outerposition') % 'XTick', .8:.2:1.6,
xlim(xl)
ylim(yl)
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS2Bii'))

%% FF PSTH
temp = [FS.FanoFactor];
temp = [temp.Sliding250ms];
hFFmat = [temp.MatchedHit];

temp = [FS.FanoFactor];
temp = [temp.Sliding250ms];
mFFmat = [temp.MatchedMiss];

lw=4;
yl=[1 2];
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
  shadedErrorBar(FS(1).ISI.param.Sliding250msTimeline, nanmean(mFFmat(:,u2plot),2), ...
    nanstd(mFFmat(:,u2plot),0,2)/sqrt(numel(u2plot)), {'r', 'LineWidth', lw},1)
  shadedErrorBar(FS(1).ISI.param.Sliding250msTimeline, nanmean(hFFmat(:,u2plot),2), ...
    nanstd(hFFmat(:,u2plot),0,2)/sqrt(numel(u2plot)), {'b', 'LineWidth', lw},1)
  set(gca, 'XTick', [-250 -125 0 100 500], 'FontSize', fs)
  xlim([-250 100])
  if spi==1
    ylabel('Fano Factor', 'FontSize', fs)
    text(100,yl(1)+.125*range(yl), 'Hit', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'blue', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    text(100,yl(1), 'Miss', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
  end
  if spi==2
    xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
  end
  ylim(yl)
  % title(' ', 'FontSize', fs)
  title(strcat(utitle,' (N=',num2str(numel(u2plot)), ')'), 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS2C'))

%% pre post stats
ff = [FS.FanoFactor];
tci=isitcind.pre250ms;

figure('Position',[100 100 900 180])
for ii=1:3
  switch ii
    case 1
      % % FF stick iqr distribution hit miss
      yl=[.5 2];
      tctitle='Pre vs Post 100ms';
      temp=[ff.Pre100ms];
      hp = [temp.Maxstim];
      temp=[ff.Post100ms];
      mp = [temp.Maxstim];
      hcol=[0 0 0];
      mcol=[1 0 1];
    case 2
      % % FF stick iqr distribution hit miss
      yl=[2.5 10];
      tctitle='Pre250ms';
      temp=[ff.Pre250ms];
      hp = [temp.MatchedHit];
      mp = [temp.MatchedMiss];
      hcol=[0 0 1];
      mcol=[1 0 0];
    case 3
      % % FF stick iqr distribution hit miss
      yl=[.5 2];
      tctitle='Post100ms';
      temp=[ff.Post100ms];
      hp = [temp.MatchedHit];
      mp = [temp.MatchedMiss];
      hcol=[0 0 1];
      mcol=[1 0 0];
  end
  
  subplot(1,3,ii)
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
    plot(nanmedian(hp(u2plot)), 2.5*spi-1, 's', 'Color', hcol, 'MarkerFaceColor', hcol, 'MarkerSize', 10)
    plot([prctile(hp(u2plot),25) prctile(hp(u2plot),75)], 2.5*[spi spi]-1, 'Color', hcol, 'LineWidth', 2.5)
    plot(nanmedian(mp(u2plot)), 2.5*spi, 's', 'Color', mcol, 'MarkerFaceColor', mcol, 'MarkerSize', 10)
    plot([prctile(mp(u2plot),25) prctile(mp(u2plot),75)], 2.5*[spi spi], 'Color', mcol, 'LineWidth', 2.5)
    pval=signrank(hp(u2plot), mp(u2plot));
    %     text(yl(2), 2.5*spi-0.5,sprintf('p=%.3f', pval*3))
    if pval<0.05/3
      plot(yl(2), 2.5*spi-0.5,'k*', 'LineWidth', 1)
    end
    ut=[ut utitle];
  end
  if ii==1
    set(gca, 'YTick', 2.5*[1:3]-0.5, 'YTickLabel', ut, 'YDir', 'reverse', 'FontSize', fs)
  else
    set(gca, 'YTick', 2.5*[1:3]-0.5, 'YTickLabel', [], 'YDir', 'reverse', 'FontSize', fs)
  end
  xlabel(tctitle, 'FontSize', fs)
  if ii==2
    title('Fano Factor', 'FontSize', fs)
  else
    title(' ', 'FontSize', fs)
  end
  xlim(yl)
  %   xlim([yl(1) yl(2)+.1*(yl(2)-yl(1))])
  ylim([1 8])
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS2D'))

%% prepre isi 2d histogram iid subtract
cb=false;
optlogscale=0;
if optlogscale==1
  xisi=log10(1000./FS(1).ISI.param.NormProbtimeline);
  ft=[8 15 30 60 120];
  xl=log10([ft(1) ft(end)]);
elseif optlogscale==0
  xisi=FS(1).ISI.param.NormProbtimeline;
  xl=[3 100];
end
yl=xl;

temp = [FS.ISI];
temp = [temp.ConsecISI];
cosecisi2dhist = cat(3, temp.hist2D);
cosecisi2dhist_iid = cat(3, temp.hist2D_ifRenewal);

figure('Position',[0 0 900 900])
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
  subplot(3,3,spi)
  hold all
  h=pcolor(xisi,xisi,squeeze(mean(cosecisi2dhist(:,:,u2plot),3)));
  set(h, 'EdgeColor', 'none')
  if cb
    colorbar
  end
  caxis([0 0.0002])
  xlim(xl)
  ylim(yl)
  if optlogscale==1
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', yt, 'FontSize', fs)
    if spi==2
      xlabel('log ISI^-^1 (s^-^1): Before', 'FontSize', fs)
    end
    if spi==1
      ylabel('log ISI^-^1 (s^-^1): After', 'FontSize', fs)
    end
  else
    set(gca, 'XTick', 0:25:100, 'YTick', 0:25:100, 'FontSize', fs)
    if spi==2
      xlabel('ISI (ms): Before', 'FontSize', fs)
    end
    if spi==1
      ylabel('ISI (ms): After', 'FontSize', fs)
    end
  end
  if spi==1
    title(strcat(utitle, ' (N=', num2str(numel(u2plot)), '): All Trials'), 'FontSize', fs)
  else
    title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
  end
  
  subplot(3,3,spi+3)
  hold all
  h=pcolor(xisi,xisi,squeeze(mean(cosecisi2dhist_iid(:,:,u2plot),3)));
  set(h, 'EdgeColor', 'none')
  if cb
    colorbar
  end
  caxis([0 0.0002])
  xlim(xl)
  ylim(yl)
  if optlogscale==1
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', log10(ft), 'YTickLabel', ft, 'FontSize', fs)
    if spi==2
      xlabel('log ISI^-^1 (s^-^1): Before', 'FontSize', fs)
    end
    if spi==1
      ylabel('log ISI^-^1 (s^-^1): After', 'FontSize', fs)
    end
  else
    set(gca, 'XTick', 0:25:100, 'YTick', 0:25:100, 'FontSize', fs)
    if spi==2
      xlabel('ISI (ms): Before', 'FontSize', fs)
    end
    if spi==1
      ylabel('ISI (ms): After', 'FontSize', fs)
    end
  end
  if spi==1
    title('If Renewal', 'FontSize', fs)
  end
  
  subplot(3,3,spi+6)
  hold all
  h=pcolor(xisi,xisi,squeeze(mean(cosecisi2dhist(:,:,u2plot)-cosecisi2dhist_iid(:,:,u2plot),3)));
  set(h, 'EdgeColor', 'none')
  if cb
    colorbar
  end
  caxis([-0.00005 0.00005])
  xlim(xl)
  ylim(yl)
  if optlogscale==1
    set(gca, 'XTick', log10(ft), 'XTickLabel', ft, 'YTick', log10(ft), 'YTickLabel', ft, 'FontSize', fs)
    if spi==2
      xlabel('log ISI^-^1 (s^-^1): Before', 'FontSize', fs)
    end
    if spi==1
      ylabel('log ISI^-^1 (s^-^1): After', 'FontSize', fs)
    end
  else
    set(gca, 'XTick', 0:25:100, 'YTick', 0:25:100, 'FontSize', fs)
    if spi==2
      xlabel('ISI (ms): Before', 'FontSize', fs)
    end
    if spi==1
      ylabel('ISI (ms): After', 'FontSize', fs)
    end
  end
  if spi==1
    title(strcat('(Data) - (If Renewal)'), 'FontSize', fs)
  end
end
% colormap jet
if cb
  print('-painters', '-dpdf', strcat(figPath, 'figS2E_cb'))
else
  print('-painters', '-dpdf', strcat(figPath, 'figS2E'))
end

%% percent change between consecutive ISIs: compare with 1000X shuffle distribution
temp = [FS.ISI];
temp = [temp.ConsecISI];
cosecisipdf = cat(1, temp.pcm_pdf);
cosecisipdfshufmean = cat(1, temp.pcm_pdf_shufflemean);
cosecisivar = cat(1, temp.pcm_var);
cosecisivarshufci = cat(1, temp.pcm_var_shuffleCI);

xl=[-2 2]; yl=[0 0.85];
figure('Position', [0 0 900 300])
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
  % plot(-2.995:.01:2.995, smooth(mean(pdfISIpcm_shuffled_mega(u2plot,:),1),5), 'k-','LineWidth', 1)
  shadedErrorBar(-2.995:.01:2.995, smooth(mean(cosecisipdfshufmean(u2plot,:),1),5), ...
    smooth(std(cosecisipdfshufmean(u2plot,:),0,1)/sqrt(numel(u2plot)),5), {'Color', 'k', 'LineWidth', 2}, 1)
  shadedErrorBar(-2.995:.01:2.995, smooth(mean(cosecisipdf(u2plot,:),1),5), ...
    smooth(std(cosecisipdf(u2plot,:),0,1)/sqrt(numel(u2plot)),5), {'Color', uc, 'LineWidth', 2}, 1)
  xlim(xl)
  ylim(yl)
  
  % shadedErrorBar(-2.995:.01:2.995, smooth(mean(pdfISIpcm_mega(u2plot,:)-pdfISIpcm_shuffled_mega(u2plot,:),1),5), smooth(std(pdfISIpcm_mega(u2plot,:),0,1)/sqrt(numel(u2plot)),5), {'Color', uc, 'LineWidth', 2}, 1)
  % ylim([-0.1 0.15])
  
  signarrow=cosecisivar(u2plot)<cosecisivarshufci(u2plot,1);
  text(xl(1), yl(2), strcat(sprintf(' %d/%d ', nnz(signarrow), numel(u2plot)), {' '}, utitle), ...
    'Color', uc, 'FontSize',fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
  text(xl(1), yl(2)-.1*range(yl), ' Observed', 'Color', uc, 'FontSize',14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
  text(xl(1), yl(2)-.2*range(yl), ' Shuffled', 'Color', 'k', 'FontSize',14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
  
  set(gca, 'FontSize',fs)
  % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize',fs)
  title(' ', 'FontSize',fs)
  if spi==1
    ylabel('Probability', 'FontSize',fs)
  elseif spi==2
    xlabel('log Percent Change in Consecutive ISIs', 'FontSize',fs)
  end
end
print('-painters', '-dpdf', strcat(figPath, 'figS2F'))
