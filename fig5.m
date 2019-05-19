sumind=struct();
for ii=1:numel(Session(1).TrialSummary_ColumnLabel)
  sumind.(Session(1).TrialSummary_ColumnLabel{ii})=ii;
end

%% psychometric curves
temp=[Session.PsychometricHR];
figure('Position', [0 0 250 300])
hold all
plot(Session(1).PsychometricStimAmp, 100*temp, 'Color', [0 0 0 0.1])
plot(Session(1).PsychometricStimAmp, 100*mean(temp,2), 'k-', 'LineWidth', 3)
xlim([0 3.6])
ylim([0 100])
set(gca, 'XTick', [0 3.6], 'XTickLabel', {'Catch', 'Max'}, 'YTick', 0:50:100, 'FontSize', fs)
xlabel('Stimulus Amplitude', 'FontSize', fs)
ylabel('Hit Rate (%)', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig5A'))

%% RT stick iqr distribution across trial types
cumalltrials=vertcat(Session.TrialSummary);

yl=[175 475];
temprt=[];
templab=[];
figure('Position',[200 200 150 300])
hold all
for ii=1:3
  switch ii
    case 0
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)==0;
    case 1
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)>0 & cumalltrials(:,sumind.stimampNormThresh)<=1;
    case 2
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)>1 & cumalltrials(:,sumind.stimbin)<16;
    case 3
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)==16;
  end
  errorbar(ii, nanmedian(1000*cumalltrials(trials2plot,sumind.RT)), ...
    prctile(1000*cumalltrials(trials2plot,sumind.RT),50)-prctile(1000*cumalltrials(trials2plot,sumind.RT),25), ...
    prctile(1000*cumalltrials(trials2plot,sumind.RT),75)-prctile(1000*cumalltrials(trials2plot,sumind.RT),50), ...
    's-', 'MarkerSize', 10, 'CapSize', 0, 'Color', ii/3*[1 0 1], 'MarkerFaceColor', ii/3*[1 0 1], 'LineWidth', 2.5)
  temprt=[temprt; cumalltrials(trials2plot,sumind.RT)];
  templab=[templab; ii*ones(nnz(trials2plot),1)];
end

[p,tbl,stats] = kruskalwallis(temprt, templab, 'off');
disp(p)
% text(.5, xl(1),sprintf('p=%.3f', p), 'FontSize', 16, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
% c=multcompare(stats,'CType','dunn-sidak','Display','off');
c=multcompare(stats,'CType','bonferroni','Display','off');
sigp=find(c(:,6)<0.05);
for ci=1:numel(sigp)
  plot(c(sigp(ci),1:2), yl(2)-range(yl)*0.05*[ci ci], 'k-', 'LineWidth', .5)
  plot(mean(c(sigp(ci),1:2)), yl(2)+range(yl)*0.02-range(yl)*0.05*[ci ci], 'k*', 'MarkerSize', 5, 'LineWidth', .5)
end

ylim(yl)
xlim([0.5 3.5])
ylabel('Reaction Time (ms)')
set(gca, 'XTick', 1:3, 'XTickLabel', {'Sub', 'Supra', 'Max'}, 'XTickLabelRotation', 45, 'YTickLabelRotation', 90, 'YTick', 0:100:700, 'FontSize', 14)
% set(gca, 'XTick', [],'FontWeight', 'bold', 'FontSize', 14)

trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1;
sprintf('%.0f %.0f %.0f', prctile(1000*cumalltrials(trials2plot,sumind.RT),50), prctile(1000*cumalltrials(trials2plot,sumind.RT),25), prctile(1000*cumalltrials(trials2plot,sumind.RT),75))
print('-painters', '-dpdf', strcat(figPath, 'fig5Bi'))

%% RT CDF
RTmedianiqr=zeros(4,3);

figure('Position', [0 0 300 300])
hold all
for ii=3:-1:1
  switch ii
    case 0
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)==0;
    case 1
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)>0 & cumalltrials(:,sumind.stimampNormThresh)<=1;
    case 2
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)>1 & cumalltrials(:,sumind.stimbin)<16;
    case 3
      trials2plot=cumalltrials(:,sumind.BehaviorOutcome)==1 & cumalltrials(:,sumind.AcceptCleanTrials)==1 & ...
        cumalltrials(:,sumind.stimbin)==16;
  end
  [f,x]=ecdf(1000*cumalltrials(trials2plot,sumind.RT));
  plot(x,f, '.-', 'Color', ii/3*[1 0 1], 'LineWidth',1.5)
  
  [fm, fi]=min(abs(f-0.5));
  RTmedianiqr(ii,1)=x(fi);
  [fm, fi]=min(abs(f-0.25));
  RTmedianiqr(ii,2)=x(fi);
  [fm, fi]=min(abs(f-0.75));
  RTmedianiqr(ii,3)=x(fi);
  
end
% plot([100 100], [0 0.47], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
plot([500 500], [0 0.46], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
plot([700 700], [0 0.67], 'Color', [0 0 0 0.1], 'LineWidth', 1.5)
text(0,0.05,'Stimulus Start', 'FontAngle', 'italic', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Rotation', 90,'FontSize', 14)
text(500,0.05,'Stimulus End', 'FontAngle', 'italic', 'VerticalAlignment', 'top', 'Rotation', 90,'FontSize', 14)
text(700,0.05,'Reward Window End', 'FontAngle', 'italic', 'VerticalAlignment', 'top', 'Rotation', 90,'FontSize', 14)
legend('boxoff')
xlim([0 700])
set(gca, 'XTick', [0 100 250 500 700], 'FontSize', fs)
xlabel('Reaction Time (ms)', 'FontSize', fs)
ylabel('Probability', 'FontSize', fs)
title('CDF', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig5Bii'))

%% hit miss histogram matching
load(strcat(dataPath, 'HSOPDR1_170301_HitMissMatching.mat'))
fs=16;

yl=[0 53];
[hitbinhist, he]=histcounts(aligned.stimbin(dprime.hitMatchCandidate), 0.5:15.5);
hitbinhist=[hitbinhist hitbinhist(end)];
[missbinhist, me]=histcounts(aligned.stimbin(dprime.missMatchCandidate), 0.5:15.5);
missbinhist=[missbinhist missbinhist(end)];
[matchbinhist, matche]=histcounts(aligned.stimbin(dprime.matchedHit), 0.5:15.5);
matchedbinhist=histcounts(aligned.stimbin(dprime.matchedHit), 0.5:16.5);
figure('Position', [0 0 300 300])
hold all
histogram(aligned.stimbin(dprime.matchedHit), 0.5:16.5, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 0.5, 'EdgeColor','none');
stairs(me, missbinhist, 'r', 'LineWidth',4)
stairs(he, hitbinhist, 'b', 'LineWidth',4)
% bar(1:16,hitbinhist,'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor','none')
% bar(1:16,missbinhist,'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor','none')
% bar(1:16,matchedbinhist,1, 'FaceColor','none', 'EdgeColor', 'k', 'LineWidth', 2.5)
text(1,48, 'Miss', 'Color', 'r', 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom')
text(15,48, 'Hit', 'Color', 'b', 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
xlim([.5 15.5])
ylim(yl)
set(gca,'XTick',[1 8 15], 'FontSize', fs)
xlabel('Submax Stimulus Bin', 'FontSize', fs)
ylabel('Trial Count', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig5C'))

%% LFP TFR
temp = [Session.LFPPower];
temp = [temp.Sliding250ms];
htfr = cat(3, temp.MatchedHit);
mtfr = cat(3, temp.MatchedMiss);

ft=[8 15 30 60 120];
figure('Position', [0 0 400 300])
hold all
lfphm=log10(htfr)-log10(mtfr);
h=pcolor(Session(1).LFPPower.param.Sliding250msTimeline,log10(Session(1).LFPPower.param.frequency.Sliding250ms), squeeze(mean(lfphm,3))');
set(h, 'EdgeColor', 'none');
xlabel('Time from Sensory Onset (ms)', 'FontSize', fs)
ylabel('log Frequency (Hz)', 'FontSize', fs)
% title('LFP Power', 'FontSize',fs)
title(strcat('LFP Power (N=', num2str(size(lfphm,3)), ')'), 'FontSize',fs)
set(gca, 'XTick', [-250:125:0 100], 'YTick', log10(ft), 'YTickLabel', ft, 'FontSize', fs)
xlim([-250 100])
ylim([log10(ft(1)) log10(ft(end))])
caxis([-.1 .1])
colorbar
print('-painters', '-dpdf', strcat(figPath, 'fig5D'))

%% LFP quintile hit rate
LFPbands = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};
temp = [Session.LFPPower];
temp = [temp.Pre250ms];
temp = [temp.QuintileHR];

yl=[25 65];
figure('Position', [0 100 300 300])
hold all
for whichband=2:length(LFPbands)
  lfpmat=[temp.(LFPbands{whichband})];
  h=errorbar(1:5, 100*mean(lfpmat,2), 100*std(lfpmat,0,2)/sqrt(size(lfpmat,2)), ...
    'LineWidth', 1.5, 'Color', 'k', 'CapSize', 0, 'MarkerSize', 8); 
  switch LFPbands{whichband}
    case 'Beta'
      h.Marker='^';
      h.LineStyle='--';
    case 'Gamma'
      h.Marker='o';
      h.LineStyle=':';
    case 'HighGamma'
      h.Marker='d';
      h.LineStyle='-';
  end
end
legend({'15-30Hz', '30-55Hz', '66-120Hz'}, 'Location', 'SouthWest')
legend('boxoff')
set(gca, 'XTick', [1 5], 'XTickLabel', {'Low', 'High'}, 'FontSize', fs)
xlim([1 5])
ylim(yl)
xlabel('Pre 250ms LFP Power Bin', 'FontSize', fs)
ylabel('Hit Rate (%)', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig5E_pre'))

%% LFP quintile hit rate
LFPbands = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};
temp = [Session.LFPPower];
temp = [temp.Post100ms];
temp = [temp.QuintileHR];

yl=[30 60];
figure('Position', [0 100 300 300])
hold all
for whichband=2:length(LFPbands)
  lfpmat=[temp.(LFPbands{whichband})];
  h=errorbar(1:5, 100*mean(lfpmat,2), 100*std(lfpmat,0,2)/sqrt(size(lfpmat,2)), ...
    'LineWidth', 1.5, 'Color', 'k', 'CapSize', 0, 'MarkerSize', 8);
  switch LFPbands{whichband}
    case 'Beta'
      h.Marker='^';
      h.LineStyle='--';
    case 'Gamma'
      h.Marker='o';
      h.LineStyle=':';
    case 'HighGamma'
      h.Marker='d';
      h.LineStyle='-';
  end
end
legend({'20-30Hz', '30-55Hz', '66-120Hz'}, 'Location', 'SouthWest')
legend('boxoff')
set(gca, 'XTick', [1 5], 'XTickLabel', {'Low', 'High'}, 'FontSize', fs)
xlim([1 5])
ylim(yl)
xlabel('Post 100ms LFP Power Bin', 'FontSize', fs)
ylabel('Hit Rate (%)', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig5E_post'))

%% PSTH
gausssigma=3;
gausswindow=21;
tempfil= exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter=tempfil/sum(tempfil);

htt='MatchedHit';
hcol=[0 0 1];
hlab='Hit';
mtt='MatchedMiss';
mcol=[1 0 0];
mlab='Miss';

temp = [FS.PSTH];
hpsth = [temp.MatchedHit];
mpsth = [temp.MatchedMiss];

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
print('-painters', '-dpdf', strcat(figPath, 'fig5F'))
