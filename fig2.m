%% BIC and AIC post K means clustering. calculate optimal number of clusters, separability score
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

isimat=[isipeakvec; cvvec]';

idxitermat=struct();
clustercenters=struct();
AICdef=struct();
for spi=1:9
  switch spi
    case 1
      u2plot=1:numel(FS);
      utitle='FS';
    case 2
      u2plot=sFS;
      utitle='sFS';
    case 3
      u2plot=nsFS;
      utitle='nsFS';
    case 4
      u2plot=find(fssrvec==1);
      utitle='srFS';
    case 5
      u2plot=find(fssrvec==0);
      utitle='nsrFS';
    case 6
      u2plot=find(fsspmwwvec==1);
      utitle='sFS_MWW';
    case 7
      u2plot=find(fsspmwwvec==0);
      utitle='nsFS_MWW';
    case 8
      u2plot=find(fsslpval<0.05);
      utitle='slFS';
    case 9
      u2plot=find(fsslpval>=0.05);
      utitle='nslFS';
  end
  normisimat=(isimat(u2plot,:)-[0 .5])./[50 2];
  
  niters=30;
  maxclusters=15;
  n=size(normisimat,1);
  dim=size(normisimat,2);
  idxitermat.(utitle)=zeros(n,maxclusters);
  clustercenters.(utitle)=NaN(dim,maxclusters,maxclusters);
  AICdef.(utitle)=zeros(maxclusters,1); % Akaike Information Criterion
  for kclusters=1:maxclusters  % number of clusters
    [tempidx,C,sumd]=kmeans(normisimat,kclusters,'MaxIter',1000, 'Replicates', niters);%, 'Display','final');  % Matlab command for k-mean clustering
    idxitermat.(utitle)(:,kclusters)=tempidx;
    clustercenters.(utitle)(:,1:kclusters,kclusters)=C';
    RSS=sum(sumd);
    AICdef.(utitle)(kclusters)=RSS+2*kclusters*dim;
%     disp(size(C))
  end
end


%% mode ISI histogram
temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
temp=[temp.ISIpeak];

figure('Position',[0 0 300 300])
hold all
hc=histcounts(temp(sFS),0:3:250);
stairs(0:3:250-3, hc, 'Color', [0.5 0.5 0.5], 'LineWidth',4)
histogram(temp(grnsFS),0:3:250, 'EdgeColor', 'none', 'FaceColor', [0 0.5 0.5], 'FaceAlpha', 0.5);
histogram(temp(plnsFS),0:3:250, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0], 'FaceAlpha', 0.5);
yl=ylim;
yl=[0 35];
% plot([1000/30 1000/30], yl, 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
% plot([1000/55 1000/55], yl, 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
ylim(yl)
xlim([3 40])
set(gca,'XTick', 0:10:150, 'FontSize', fs)
xlabel('ISI Peak (ms)', 'FontSize', fs)
ylabel('Count', 'FontSize', fs)
text(40.5,yl(2), 'Sensory FS', 'Color', [0.5 0.5 0.5], 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
text(40.5,yl(2)-.1*(yl(2)-yl(1)), 'Gamma Regular nsFS', 'Color', [0 0.5 0.5], 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
text(40.5,yl(2)-.2*(yl(2)-yl(1)), 'Poisson-like nsFS', 'Color', [0.5 0.5 0], 'FontSize', 14, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
print('-painters', '-dpdf', strcat(figPath, 'fig2A'))

%% ISI PDF FS subtype averaged
temp = [FS.ISI];
temp=[temp.Pre1s];
temp=[temp.AllTrials];
temp=[temp.ISIpdf];

figure('Position', [400 100 300 300])
hold all
for spi=1:3
  switch spi
    case 1
      u2plot=sFS;
      utitle='Sensory FS';
      uc=[.5 .5 .5];
      lw=4;
    case 2
      u2plot=grnsFS;
      utitle='Gamma Non-Sensory FS';
      uc=[0 .5 .5];
      lw=4;
    case 3
      u2plot=plnsFS;
      utitle='Poisson-Like Non-Sensory FS';
      uc=[0.5 0.5 0];
      lw=4;
  end
  lw=2.5;
  plot(FS(1).ISI.param.NormProbtimeline, mean(temp(:,u2plot),2),'Color',uc, 'LineWidth', lw)
end
xlim([0 100])
yl=[0 0.02];
ylim(yl)
set(gca, 'XTick', [0:25:100], 'YTick', 0:0.01:1, 'FontSize', fs)
xlabel('ISI (ms)', 'FontSize', fs)
ylabel('Probability', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig2B'))


%% CV stick IQR horizontal + histogram
temp = [FS.ISI];
temp = [temp.Pre1s];
temp = [temp.AllTrials];
temp = [temp.CV];

he=0:0.05:4.5;
xl=[0.8 1.5];
% xlab='CV';
xlab='Coefficient of Variation';

xlabs=[];
tempcorr=[];
templab=[];

figure('Position', [300 200 270 300])
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
  tempcorr = [tempcorr temp(u2plot)];
  templab = [templab spi*ones(1,numel(u2plot))];
  
  plot(median(temp(u2plot)), spi, 's', 'MarkerEdgeColor', uc, 'MarkerFaceColor', uc, 'MarkerSize', 10)
  plot([prctile(temp(u2plot),25) prctile(temp(u2plot),75)], [spi spi], '-', 'Color', uc, 'LineWidth', 2.5)
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
set(gca, 'XTick', [], 'YTick', [], 'FontSize', fs)
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
    hc=histcounts(temp(u2plot), he);
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(temp(u2plot), he, 'EdgeColor', 'none', 'FaceColor', uc, 'FaceAlpha', .5)
  end
end
set(gca, 'XTick', .8:.2:1.6, 'FontSize', fs, 'ActivePositionProperty', 'outerposition')
xlim(xl)
ylim([0 35])
xlabel(xlab, 'FontSize', fs)
% ylabel('Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig2C'))

%% AIC
yl=[5 22];
ut=[];
figure('Position', [0 0 300 300])
hold all
for spi=1:2
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
      uc=[.5 .5 .5];
    case 2
      u2plot=nsFS;
      utitle='nsFS';
      uc=[0 .5 0];
  end
  
  plot(1:maxclusters, AICdef.(utitle), 'o-', 'MarkerFaceColor', uc, 'Color', uc, 'LineWidth', 2)
  
  ut=[ut {utitle}];
  text(5, 10-.1*range(yl)*(spi-1), utitle, 'Color', uc, 'FontSize', fs, 'FontWeight', 'bold','HorizontalAlignment', 'right','VerticalAlignment', 'top')
end
xlim([1 5])
ylim(yl)
set(gca, 'XTick', 1:10, 'FontSize', fs)
xlabel('Number of Clusters', 'FontSize', fs)
ylabel('Akaike Information Criterion ', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig2D'))

%% K means clustering
for spi=1:2
  switch spi
    case 1
      u2plot=sFS;
      utitle='sFS';
      tempucmat=[.5 .5 .5];
    case 2
      u2plot=nsFS;
      utitle='nsFS';
      tempucmat=[0 .5 .5; .5 .5 0];
  end
  
  [kclusternum,iternum]=find(AICdef.(utitle)==min(AICdef.(utitle)(:)));
  tempidx=squeeze(idxitermat.(utitle)(:,kclusternum(1),iternum(1)));
  tempuid=unique(tempidx);

  figure('Position', [300*spi 100 300 300])
  hold all
  for ui= 1:numel(tempuid)
    scatter(isipeakvec(u2plot(tempidx==tempuid(ui))), cvvec(u2plot(tempidx==tempuid(ui))) ...
      ,30, 'o', 'MarkerFaceColor', tempucmat(tempuid(ui),:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .4)
  end
  set(gca, 'YTick', 0:.2:2, 'FontSize', fs)
  xlabel('ISI Peak (ms)', 'FontSize', fs)
  ylabel('CV', 'FontSize', fs)
  axis([0 50 0.8 1.6])
  title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)

  print('-painters', '-dpdf', strcat(figPath, 'fig2E_', utitle))
end

