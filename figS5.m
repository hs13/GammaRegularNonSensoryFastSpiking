valboutunits = find(~cellfun(@isempty, {FS.Locomotion}));

%%
xl=[-1500 1500];
yl=[0 35];
gausssigma=3;
gausswindow=21;
tempfil= exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter=tempfil/sum(tempfil);

templm = [FS.Locomotion];

figure('Position', [0 100 700 300])
for ii=1:2
  switch ii
    case 1
      runspikepsthtimeline = FS(valboutunits(1)).Locomotion.param.BoutStartAligned_SpikeTrainTimeline;
      runspikepsth=conv2([templm.BoutStartAligned_SpikeTrainAvg], tempfilter, 'same');
      xlab='Time from Walking Onset (ms)';
      sptit='All Walking Bouts';
    case 2
      runspikepsthtimeline = FS(valboutunits(1)).Locomotion.param.BoutEndAligned_SpikeTrainTimeline;
      runspikepsth=conv2([templm.BoutEndAligned_SpikeTrainAvg], tempfilter, 'same');
      xlab='Time from Walking Offset (ms)';
      sptit='All Walking Bouts';
  end
  toi = runspikepsthtimeline>=-1500 & runspikepsthtimeline<1500;
    
subplot(1,2,ii)
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
  u2plot=find(ismember(valboutunits,u2plot));
  shadedErrorBar(runspikepsthtimeline(toi), nanmean(runspikepsth(toi,u2plot),2), ...
    nanstd(runspikepsth(toi,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', uc, 'LineWidth', 2.5}, 1)
if ii==2
text(xl(1), yl(2)-.1*range(yl)*(spi-1), strcat({' '},utitle, ' (N=', num2str(numel(u2plot)), ')'), ...
  'Color', uc, 'Fontsize', 14, 'Fontweight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
end
end
set(gca, 'FontSize', fs)
xlim(xl)
ylim(yl)
xlabel(xlab, 'FontSize', fs)
if ii==1
ylabel('Rate (Hz)', 'FontSize', fs)
end
title(' ', 'FontSize', fs)
% title(sptit, 'FontSize', 12)
end
print('-painters', '-dpdf', strcat(figPath, 'figS5A'))

%%
figure('Position', [0 100 700 300])
for ii=1:2
  switch ii
    case 1
      temptl=FS(valboutunits(1)).Locomotion.param.BoutStartAligned_CVTimeline;
      runcvpsth=[templm.BoutStartAligned_CV];
      xlab='Time from Walking Onset (ms)';
      sptit='All Walking Bouts';
    case 2
      temptl=FS(valboutunits(1)).Locomotion.param.BoutEndAligned_CVTimeline;
      runcvpsth=[templm.BoutEndAligned_CV];
      xlab='Time from Walking Offset (ms)';
      sptit='All Walking Bouts';
  end
  subplot(1,2,ii)
  hold all
  for spi=3:-1:1
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
        uc=[0.5 0.5 0 0.5];
    end
    u2plot=find(ismember(valboutunits,u2plot));
    shadedErrorBar(temptl, mean(runcvpsth(:,u2plot),2), std(runcvpsth(:,u2plot),0,2)/sqrt(numel(u2plot)), {'Color', uc, 'LineWidth', 4}, 1)
  end
  set(gca, 'FontSize', fs)
  xlim([-1500 1500])
  % ylim([0 30])
  xlabel(xlab, 'FontSize', fs)
  if ii==1
    ylabel('CV', 'FontSize', fs)
  end
  title(' ', 'FontSize', fs)
end
print('-painters', '-dpdf', strcat(figPath, 'figS5C'))

%% run delta rate stick IQR horizontal + histogram
pltopt=2;
switch pltopt
  case 1
    runspikepsthtimeline = FS(valboutunits(1)).Locomotion.param.BoutStartAligned_SpikeTrainTimeline;
    runspikepsth=conv2([templm.BoutStartAligned_SpikeTrainAvg], tempfilter, 'same');
    tempvar=mean(runspikepsth(runspikepsthtimeline>-500 & runspikepsthtimeline<=500,:),1) - ...
      mean(runspikepsth(runspikepsthtimeline>-1500 & runspikepsthtimeline<=-500,:),1);
    he=-15:15;
    xl=[-10 15];
    xlab='\DeltaRate (Hz)';
    xtit = 'B';
  case 2
    temptl=FS(valboutunits(1)).Locomotion.param.BoutStartAligned_CVTimeline;
    runcvpsth=[templm.BoutStartAligned_CV];
    tempvar=runcvpsth(temptl==-1,:) - runcvpsth(temptl==-1001,:);
    he=-2:0.05:2;
    xl=[-.3 .3];
    xlab='\DeltaCV';
    xtit = 'D';
end

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [200 200 250 300])
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
  u2plot=find(ismember(valboutunits,u2plot));
  xlabs=[xlabs {utitle}];
%   xlabs=[xlabs strcat(utitle, {' N='}, num2str(numel(u2plot)))];
  tempcorr = [tempcorr tempvar(u2plot)];
  templab = [templab spi*ones(1,numel(u2plot))];
  
  plot(median(tempvar(u2plot)), spi, 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 10)
  plot([prctile(tempvar(u2plot),25) prctile(tempvar(u2plot),75)], [spi spi], '-', 'Color', uc, 'LineWidth', 2.5)
  
end
[p,tbl,stats] = kruskalwallis(tempcorr', templab', 'off');
% disp(p)
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(xl(2)-range(xl)*0.05*[ci ci], c(sigp(ci),1:2), 'k-', 'LineWidth', .5)
  plot(xl(2)+range(xl)*0.02-range(xl)*0.05*[ci ci], mean(c(sigp(ci),1:2)), 'k*', 'MarkerSize', 5, 'LineWidth', .5)
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
  u2plot=find(ismember(valboutunits,u2plot));
  if spi==1
    hc=histcounts(tempvar(u2plot), he);
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(tempvar(u2plot), he, 'EdgeColor', 'none', 'FaceColor', uc, 'FaceAlpha', .5)
  end
end
set(gca, 'FontSize', fs, 'ActivePositionProperty', 'outerposition')
xlim(xl)
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'figS5', xtit))

