temp=[FS.SensoryAltDef];
fssrvec=[temp.PostPre];
fsspmwwvec=[temp.SP_MWWtail];

temp=[FS.FirstSpikeLatency];
fsslpval = [temp.Maxstim_BrownForsythePval];

temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
isipeakvec=[temp.ISIpeak];

temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
cvvec=[temp.CV];

%%
AICdef21=struct();
for pltopt=1:4
  switch pltopt
    case 1
      nsfs2plot=nsFS;
      nsfsdeftitle='nsFS';
      grnsfsvarname='grnsFS';
      plnsfsvarname='plnsFS';
    case 2
      nsfs2plot=find(fssrvec==0);
      nsfsdeftitle='nsrFS';
      grnsfsvarname='grnsrFS';
      plnsfsvarname='plnsrFS';
    case 3
      nsfs2plot=find(fsspmwwvec==0);
      nsfsdeftitle='nsFS_MWW';
      grnsfsvarname='grnsFS_MWW';
      plnsfsvarname='plnsFS_MWW';
    case 4
      nsfs2plot=find(fsslpval>=0.05);
      nsfsdeftitle='nslFS';
      grnsfsvarname='grnslFS';
      plnsfsvarname='plnslFS';
  end
  tempidx=squeeze(idxitermat.(nsfsdeftitle)(:,2)); % 2 clusters, 1st iteration
  tempcc=squeeze(clustercenters.(nsfsdeftitle)(:,1:2,2));
  if tempcc(1)<tempcc(2)
    AICdef21.(plnsfsvarname)=nsfs2plot(tempidx==1);
    AICdef21.(grnsfsvarname)=nsfs2plot(tempidx==2);
  else
    AICdef21.(plnsfsvarname)=nsfs2plot(tempidx==2);
    AICdef21.(grnsfsvarname)=nsfs2plot(tempidx==1);
  end
end

isequal(grnsFS, AICdef21.grnsFS)
isequal(plnsFS, AICdef21.plnsFS)

%%
yl=[0 25];
figure('Position', [0 100 1200 330])
for pltopt=1:4
  switch pltopt
    case 1
      sfs2plot=sFS;
      nsfs2plot=nsFS;
      sfsdeftitle='sFS';
      nsfsdeftitle='nsFS';
      sfstitle='sFS';
      nsfstitle='nsFS';
      subtitle={'Max vs Catch', 'Ideal Observer'};
    case 2
      sfs2plot=find(fssrvec==1);
      nsfs2plot=find(fssrvec==0);
      sfsdeftitle='srFS';
      nsfsdeftitle='nsrFS';
      sfstitle='srFS';
      nsfstitle='nsrFS';
      subtitle={'Post vs Pre', 'Wilcoxon Signed Rank  '};
    case 3
      sfs2plot=find(fsspmwwvec==1);
      nsfs2plot=find(fsspmwwvec==0);
      sfsdeftitle='sFS_MWW';
      nsfsdeftitle='nsFS_MWW';
      sfstitle='sFS-MWW';
      nsfstitle='nsFS-MWW';
      subtitle={'Max vs Catch', '  Mann-Whitney-Wilcoxon'};
    case 4
      sfs2plot=find(fsslpval<0.05);
      nsfs2plot=find(fsslpval>=0.05);
      sfsdeftitle='slFS';
      nsfsdeftitle='nslFS';
      sfstitle='slFS';
      nsfstitle='nslFS';      
      subtitle={'1^s^t Spike Latency', 'Brown-Forsythe'};
  end
  
  subplot(1,4,pltopt)
  hold all
  for spi=1:2
    switch spi
      case 1
        u2plot=sfs2plot;
        utitle=sfstitle;
        udeftitle=sfsdeftitle;
        uc=[.5 .5 .5];
      case 2
        u2plot=nsfs2plot;
        utitle=nsfstitle;
        udeftitle=nsfsdeftitle;
        uc=[0 .5 0];
    end
    
    for ii=1:size(AICdef.(udeftitle),2)
      scatter(1:maxclusters, AICdef.(udeftitle)(:,ii) ...
        , 40, 'o', 'MarkerFaceColor', uc, 'MarkerEdgeColor', 'none')
    end
    % plot(1:maxclusters, min(AIC.(utitle),[],2), 'r-', 'LineWidth', 2)
    plot(1:maxclusters, mean(AICdef.(udeftitle),2), 'Color', uc, 'LineWidth', 2)
    %   ut=[ut {utitle}];
    if pltopt==0
      text(5, yl(1)+.2*range(yl)-.1*range(yl)*spi, strcat(num2str(numel(sfs2plot)), {' '}, utitle), ...
        'Color', uc, 'FontSize', 14, 'FontWeight', 'bold','HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
    else
      if spi==1
        text(5, yl(1)+.1*range(yl), strcat(sprintf('%d / %d ', nnz(ismember(sfs2plot, sFS)), numel(sfs2plot)), {' '}, utitle), ...
          'Color', uc, 'FontSize', 14, 'FontWeight', 'bold','HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
      elseif spi==2
        text(5, yl(1), strcat(sprintf('%d / %d ', nnz(ismember(nsfs2plot, nsFS)), numel(nsfs2plot)), {' '}, utitle), ...
          'Color', uc, 'FontSize', 14, 'FontWeight', 'bold','HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
      end
    end
  end
  xlim([1 5])
  ylim(yl)
  % plot(2:maxclusters, diff(min(AIC.(utitle),[],2)), 'r-', 'LineWidth', 2)
  % plot(2:maxclusters, diff(mean(AIC.(utitle),2)), 'c-', 'LineWidth', 2)
  set(gca, 'XTick', 1:10, 'FontSize', fs)
  xlabel('Number of Clusters', 'FontSize', fs)
  if pltopt==1
    ylabel('AIC', 'FontSize', fs)
  end
  title(subtitle, 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS1A'))

%% K means clustering
xl=[0 50]; yl=[0.8 1.6];
figure('Position', [0 100 1200 300])
for pltopt=1:4
  switch pltopt
    case 1
      sfs2plot=sFS;
      nsfs2plot=nsFS;
      sfsdeftitle='sFS';
      nsfsdeftitle='nsFS';
      sfstitle='sFS';
      nsfstitle='nsFS';
    case 2
      sfs2plot=find(fssrvec==1);
      nsfs2plot=find(fssrvec==0);
      sfsdeftitle='srFS';
      nsfsdeftitle='nsrFS';
      sfstitle='srFS';
      nsfstitle='nsrFS';
    case 3
      sfs2plot=find(fsspmwwvec==1);
      nsfs2plot=find(fsspmwwvec==0);
      sfsdeftitle='sFS_MWW';
      nsfsdeftitle='nsFS_MWW';
      sfstitle='sFS-MWW';
      nsfstitle='nsFS-MWW';
    case 4
      sfs2plot=find(fsslpval<0.05);
      nsfs2plot=find(fsslpval>=0.05);
      sfsdeftitle='slFS';
      nsfsdeftitle='nslFS';
      sfstitle='slFS';
      nsfstitle='nslFS';
  end
  
  u2plot=nsfs2plot;
  udeftitle=nsfsdeftitle;
  utitle=nsfstitle;
  
  %   [kclusternum,iternum]=find(AICdef.(utitle)==min(AICdef.(utitle)(:)));
  %   tempidx=squeeze(idxitermat.(utitle)(:,kclusternum(1),iternum(1)));
  tempidx=idxitermat.(udeftitle)(:,2); % 2 clusters, 1st iteration
  tempuid=unique(tempidx);
  tempcc=squeeze(clustercenters.(udeftitle)(:,1:2,2));
  if tempcc(1,1)<tempcc(1,2)
    tempucmat=[.5 .5 0; 0 .5 .5];
  else
    tempucmat=[0 .5 .5; .5 .5 0];
  end
  
  subplot(1,4,pltopt)
  hold all
  for ui= 1:numel(tempuid)
    scatter(isipeakvec(u2plot(tempidx==tempuid(ui))), cvvec(u2plot(tempidx==tempuid(ui))) ...
      ,30, 'o', 'MarkerFaceColor', tempucmat(tempuid(ui),:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .4)    
  end
  
  text(xl(2), yl(2)-.1*range(yl), strcat('pl', utitle, ' (N=', num2str(nnz(tempidx== (tempcc(1,1)>tempcc(1,2))+1)), ')'), ...
    'FontSize', 14, 'Color', tempucmat((tempcc(1,1)>tempcc(1,2))+1,:), 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
  text(xl(2), yl(2), strcat('gr', utitle, ' (N=', num2str(nnz(tempidx== (tempcc(1,1)<tempcc(1,2))+1)), ')'), ...
    'FontSize', 14, 'Color', tempucmat((tempcc(1,1)<tempcc(1,2))+1,:), 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
  
  set(gca, 'YTick', 0:.2:2, 'FontSize', fs)
%   if pltopt==2
    xlabel('ISI Peak (ms)', 'FontSize', fs)
%   end
  if pltopt==1
    ylabel('CV', 'FontSize', fs)
  end
  axis([xl yl])
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS1B'))

%% mode ISI / CV stick IQR horizontal + histogram
varopt=1;
switch varopt
  case 1
    tempvar=isipeakvec;
    he=0:3:250;
    xl=[3 40];
    xlab='ISI Peak (ms)';
    xtit = 'C';
  case 2
    tempvar=cvvec;
    he=0:0.05:4.5;
    xl=[0.8 1.5];
    xlab='CV';
    % xlab='Coefficient of Variation';
    xtit = 'D';
end

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [0 100 1200 300])
for pltopt=1:4
  switch pltopt
    case 1
      sfs2plot=sFS;
      nsfs2plot=nsFS;
      grnsfs2plot=AICdef21.grnsFS;
      plnsfs2plot=AICdef21.plnsFS;
      sfsdeftitle='sFS';
      nsfsdeftitle='nsFS';
      sfstitle='sFS';
      nsfstitle='nsFS';
      grnsfstitle='grnsFS';
      plnsfstitle='plnsFS';
    case 2
      sfs2plot=find(fssrvec==1);
      nsfs2plot=find(fssrvec==0);
      grnsfs2plot=AICdef21.grnsrFS;
      plnsfs2plot=AICdef21.plnsrFS;
      sfsdeftitle='srFS';
      nsfsdeftitle='nsrFS';
      sfstitle='srFS';
      nsfstitle='nsrFS';
      grnsfstitle='grnsrFS';
      plnsfstitle='plnsrFS';
    case 3
      sfs2plot=find(fsspmwwvec==1);
      nsfs2plot=find(fsspmwwvec==0);
      grnsfs2plot=AICdef21.grnsFS_MWW;
      plnsfs2plot=AICdef21.plnsFS_MWW;
      sfsdeftitle='sFS_MWW';
      nsfsdeftitle='nsFS_MWW';
      sfstitle='sFS-MWW';
      nsfstitle='nsFS-MWW';
      grnsfstitle='grnsFS-MWW';
      plnsfstitle='plnsFS-MWW';
    case 4
      sfs2plot=find(fsslpval<0.05);
      nsfs2plot=find(fsslpval>=0.05);
      grnsfs2plot=AICdef21.grnslFS;
      plnsfs2plot=AICdef21.plnslFS;
      sfsdeftitle='slFS';
      nsfsdeftitle='nslFS';
      sfstitle='slFS';
      nsfstitle='nslFS';
      grnsfstitle='grnslFS';
      plnsfstitle='plnslFS';
  end
  
subplot(5,4,pltopt)
hold all
for spi=1:3
    switch spi
      case 1
        u2plot=sfs2plot;
        utitle=sfstitle;
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsfs2plot;
        utitle=grnsfstitle;
        uc=[0 .5 .5];
      case 3
        u2plot=plnsfs2plot;
        utitle=plnsfstitle;
        uc=[.5 .5 0];
    end
  xlabs=[xlabs {utitle}];
%   xlabs=[xlabs strcat(utitle, {' N='}, num2str(numel(u2plot)))];
  tempcorr = [tempcorr tempvar(u2plot)];
  templab = [templab spi*ones(1,numel(u2plot))];
  
  plot(median(tempvar(u2plot)), spi, 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 10)
  plot([prctile(tempvar(u2plot),25) prctile(tempvar(u2plot),75)], [spi spi], '-', 'Color', uc, 'LineWidth', 2.5)
  
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

yl=[0 35];
subplot(5,4,4*[1:4]+pltopt)
hold all
for spi=1:3
    switch spi
      case 1
        u2plot=sfs2plot;
        utitle=sfstitle;
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsfs2plot;
        utitle=grnsfstitle;
        uc=[0 .5 .5];
      case 3
        u2plot=plnsfs2plot;
        utitle=plnsfstitle;
        uc=[.5 .5 0];
    end
  
  if spi==1
    hc=histcounts(tempvar(u2plot), he);
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(tempvar(u2plot), he, 'EdgeColor', 'none', 'FaceColor', uc, 'FaceAlpha', .5)
  end
  if varopt==1
  text(xl(2), yl(2)-.1*range(yl)*(spi-1), utitle, ...
    'Color', uc, 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
  end
end
set(gca, 'XTick', .8:.2:1.6, 'FontSize', fs, 'ActivePositionProperty', 'outerposition')
xlim(xl)
ylim(yl)
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS1', xtit))

%% neurometric spike count FS
temp=[FS.SpikeRate];
temp=[temp.Post100ms];
thrmat=[[temp.Catch]; [temp.Sub]; [temp.Supra]; [temp.Maxstim]];

yl=[8 25];
figure('Position', [0 100 1200 300])
for pltopt=1:4
  switch pltopt
    case 1
      sfs2plot=sFS;
      nsfs2plot=nsFS;
      grnsfs2plot=AICdef21.grnsFS;
      plnsfs2plot=AICdef21.plnsFS;
      sfsdeftitle='sFS';
      nsfsdeftitle='nsFS';
      sfstitle='sFS';
      nsfstitle='nsFS';
      grnsfstitle='grnsFS';
      plnsfstitle='plnsFS';
    case 2
      sfs2plot=find(fssrvec==1);
      nsfs2plot=find(fssrvec==0);
      grnsfs2plot=AICdef21.grnsrFS;
      plnsfs2plot=AICdef21.plnsrFS;
      sfsdeftitle='srFS';
      nsfsdeftitle='nsrFS';
      sfstitle='srFS';
      nsfstitle='nsrFS';
      grnsfstitle='grnsrFS';
      plnsfstitle='plnsrFS';
    case 3
      sfs2plot=find(fsspmwwvec==1);
      nsfs2plot=find(fsspmwwvec==0);
      grnsfs2plot=AICdef21.grnsFS_MWW;
      plnsfs2plot=AICdef21.plnsFS_MWW;
      sfsdeftitle='sFS_MWW';
      nsfsdeftitle='nsFS_MWW';
      sfstitle='sFS-MWW';
      nsfstitle='nsFS-MWW';
      grnsfstitle='grnsFS-MWW';
      plnsfstitle='plnsFS-MWW';
    case 4
      sfs2plot=find(fsslpval<0.05);
      nsfs2plot=find(fsslpval>=0.05);
      grnsfs2plot=AICdef21.grnslFS;
      plnsfs2plot=AICdef21.plnslFS;
      sfsdeftitle='slFS';
      nsfsdeftitle='nslFS';
      sfstitle='slFS';
      nsfstitle='nslFS';
      grnsfstitle='grnslFS';
      plnsfstitle='plnslFS';
  end
  
  subplot(1,4,pltopt)
  hold all
  for spi=3:-1:1
    switch spi
      case 1
        u2plot=sfs2plot;
        utitle=sfstitle;
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsfs2plot;
        utitle=grnsfstitle;
        uc=[0 .5 .5];
      case 3
        u2plot=plnsfs2plot;
        utitle=plnsfstitle;
        uc=[.5 .5 0];
    end
    lw=2.5;
    plot(1:size(thrmat,1), nanmean(thrmat(:,u2plot),2), 'Color', uc, 'LineWidth', lw)
    for ii=1:size(thrmat,1)
      errorbar(ii, nanmean(thrmat(ii,u2plot),2), nanstd(thrmat(ii,u2plot),0,2)/sqrt(numel(u2plot)), ...
        'o-', 'MarkerSize', 5, 'Color', uc, 'MarkerFaceColor', uc, 'LineWidth', lw)
    end
  end
%   set(gca, 'XTick', 1:size(thrmat,1), 'XTickLabel', {'catch', 'sub', 'supra', 'maxstim'}, 'XTickLabelRotation', 45, 'FontSize', fs)%, 'ActivePositionProperty', 'outerposition')
  set(gca, 'XTick', [1 size(thrmat,1)], 'XTickLabel', {'Catch', 'Max'}, 'FontSize', fs)
  ylim(yl)
  xlim([0.5 size(thrmat,1)+0.5])
  xlabel('Stimulus Amplitude', 'FontSize', fs)
  if pltopt==1
    ylabel('Rate (Hz)', 'FontSize', fs)
  end
  % title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS1E'))

%% spike latency
temp = [FS.FirstSpikeLatency];
temp = [temp.Maxstim];
slmat = [temp.PDF];

figure('Position', [0 100 1200 300])
for pltopt=1:4
  switch pltopt
    case 1
      sfs2plot=sFS;
      nsfs2plot=nsFS;
      grnsfs2plot=AICdef21.grnsFS;
      plnsfs2plot=AICdef21.plnsFS;
      sfsdeftitle='sFS';
      nsfsdeftitle='nsFS';
      sfstitle='sFS';
      nsfstitle='nsFS';
      grnsfstitle='grnsFS';
      plnsfstitle='plnsFS';
    case 2
      sfs2plot=find(fssrvec==1);
      nsfs2plot=find(fssrvec==0);
      grnsfs2plot=AICdef21.grnsrFS;
      plnsfs2plot=AICdef21.plnsrFS;
      sfsdeftitle='srFS';
      nsfsdeftitle='nsrFS';
      sfstitle='srFS';
      nsfstitle='nsrFS';
      grnsfstitle='grnsrFS';
      plnsfstitle='plnsrFS';
    case 3
      sfs2plot=find(fsspmwwvec==1);
      nsfs2plot=find(fsspmwwvec==0);
      grnsfs2plot=AICdef21.grnsFS_MWW;
      plnsfs2plot=AICdef21.plnsFS_MWW;
      sfsdeftitle='sFS_MWW';
      nsfsdeftitle='nsFS_MWW';
      sfstitle='sFS-MWW';
      nsfstitle='nsFS-MWW';
      grnsfstitle='grnsFS-MWW';
      plnsfstitle='plnsFS-MWW';
    case 4
      sfs2plot=find(fsslpval<0.05);
      nsfs2plot=find(fsslpval>=0.05);
      grnsfs2plot=AICdef21.grnslFS;
      plnsfs2plot=AICdef21.plnslFS;
      sfsdeftitle='slFS';
      nsfsdeftitle='nslFS';
      sfstitle='slFS';
      nsfstitle='nslFS';
      grnsfstitle='grnslFS';
      plnsfstitle='plnslFS';
  end
  
  subplot(1,4,pltopt)
  hold all
  for spi=1:3
    switch spi
      case 1
        u2plot=sfs2plot;
        utitle=sfstitle;
        uc=[.5 .5 .5];
      case 2
        u2plot=grnsfs2plot;
        utitle=grnsfstitle;
        uc=[0 .5 .5];
      case 3
        u2plot=plnsfs2plot;
        utitle=plnsfstitle;
        uc=[.5 .5 0];
    end
    hold all
    plot(FS(1).FirstSpikeLatency.param.binnedSpikeLatencyNormProbtimeline, mean(slmat(:,u2plot),2), ...
      'Color', uc, 'LineWidth', 2.5)
    xlim([0 50])
    yl=[0 0.028];
    ylim(yl)
    set(gca, 'XTick', 0:25:100, 'YTick', [0 0.025],'FontSize', fs)%, 'ActivePositionProperty', 'outerposition')
      xlabel('1^s^t Spike Latency (ms)', 'FontSize', fs)
    if pltopt==1
      ylabel('Probability', 'FontSize', fs)
    end
  end
end
orient(gcf, 'landscape')
print('-painters', '-dpdf', strcat(figPath, 'figS1F'))
