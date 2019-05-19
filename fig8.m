load(strcat(dataPath, 'exampleRS.mat'))
t0ind=find(exRS.TrialSpikeTrainTimeline==0);
s2a=find(strcmp({Session.ID}, exRS.session));
trialSummary = Session(s2a).TrialSummary;

imh=find(Session(s2a).TrialSummary(:,sumind.MatchedTrials)==1 & Session(s2a).TrialSummary(:,sumind.BehaviorOutcome)==1);
[~,ih]=sort(Session(s2a).TrialSummary(imh,sumind.stimamp));
hsup=find(Session(s2a).TrialSummary(imh(ih),sumind.stimampNormThresh)>=1, 1, 'first');
[hspkt,htri]=find(exRS.TrialSpikeTrain(:,imh(ih))==1);
[hspkt1i,h1sttri]=find(cumsum(exRS.TrialSpikeTrain(t0ind:end,imh(ih)),1)==1 & exRS.TrialSpikeTrain(t0ind:end,imh(ih))==1);


imm=find(Session(s2a).TrialSummary(:,sumind.MatchedTrials)==1 & Session(s2a).TrialSummary(:,sumind.BehaviorOutcome)==0);
[~,im]=sort(Session(s2a).TrialSummary(imm,sumind.stimamp));
msup=find(Session(s2a).TrialSummary(imm(im),sumind.stimampNormThresh)>=1, 1, 'first');
[mspkt,mtri]=find(exRS.TrialSpikeTrain(:,imm(im))==1);
[mspkt1i,m1sttri]=find(cumsum(exRS.TrialSpikeTrain(t0ind:end,imm(im)),1)==1 & exRS.TrialSpikeTrain(t0ind:end,imm(im))==1);

% %%
fs=15; ms=50;
xl=[-50 50];
xt=-50:50:50;

% xl=[-250 100];
% xt=[-250 -125 0 100 500];

figure('Position',[300 300 400 300])
% figure('Position',[300 300 600 320])
subplot(1,2,1)
hold all
scatter(exRS.TrialSpikeTrainTimeline(hspkt),htri, ms, 'k.')
scatter(exRS.TrialSpikeTrainTimeline(t0ind+hspkt1i-1), h1sttri, ms, 'b.')
text(xl(1)-.1*range(xl), (hsup-0.5+numel(imh))/2, 'Suprathresh', 'Rotation', 90, 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
text(xl(1)-.1*range(xl), (hsup-0.5)/2, 'Subthresh', 'Rotation', 90, 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
set(gca, 'XTick', xt, 'YTick', hsup-0.5, 'YTickLabel', [], 'FontSize', fs)
xlim(xl)
ylim([0 numel(imh)])
xlabel(' ', 'FontSize', fs)
title('Hit', 'FontSize', fs)

subplot(1,2,2)
hold all
scatter(exRS.TrialSpikeTrainTimeline(mspkt),mtri, ms, 'k.')
scatter(exRS.TrialSpikeTrainTimeline(t0ind+mspkt1i-1), m1sttri, ms, 'r.')
text(xl(1)-.1*range(xl), (msup-0.5+numel(imm))/2, 'Suprathresh', 'Rotation', 90, 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
text(xl(1)-.1*range(xl), (msup-0.5)/2, 'Subthresh', 'Rotation', 90, 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
set(gca, 'XTick', xt, 'YTick', hsup-0.5, 'YTickLabel', [], 'FontSize', fs)
xlim(xl)
ylim([0 numel(imm)])
xlabel(' ', 'FontSize', fs)
title('Miss', 'FontSize', fs)

annotation('textbox', [0 0 1 0.08], 'String', 'Time from Sensory Onset (ms)', 'HorizontalAlignment', 'center', ...
  'VerticalAlignment', 'middle', 'FontWeight', 'normal', 'FontSize', fs, 'EdgeColor', 'none')
print('-painters', '-dpdf', strcat(figPath, 'fig8A'))

%% example cell latency histogram
fs=16;
yl=[0 0.05];
thrttalias = {'SubthreshHit', 'SubthreshMiss', 'SuprathreshHit', 'SuprathreshMiss'};

figure('Position', [100 100 300 300])
hold all
for typi=1:numel(thrttalias)
  switch thrttalias{typi}
    case 'SubthreshHit'
      ttcol = [0.4 0.4 1];
      lw = .5;
    case 'SubthreshMiss'
      ttcol = [1 0.4 0.4];
      lw = .5;
    case 'SuprathreshHit'
      ttcol = [0 0 1];
      lw = 2.5;
    case 'SuprathreshMiss'
      ttcol = [1 0 0];
      lw = 2.5;
  end
plot(exRS.FirstSpikeLatency.param.binnedSpikeLatencyNormProbtimeline, exRS.FirstSpikeLatency.(thrttalias{typi}).PDF, 'Color', ttcol, 'LineWidth', lw)
scatter(exRS.FirstSpikeLatency.(thrttalias{typi}).Timing, exRS.FirstSpikeLatency.(thrttalias{typi}).Reliability, 100, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1)
end

text(0, yl(1)+.1*range(yl), ' Hit', 'FontSize', fs, 'Color', 'b', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom')
text(0, yl(1), ' Miss', 'FontSize', fs, 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom')

ylim(yl)
xlim([0 50])
set(gca, 'XTick', -50:50:50, 'YTick', 0:0.01:1, 'FontSize', fs)
xlabel('1^s^t Spike Latency (ms)', 'FontSize', fs)
ylabel('Probability', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig8B'))

%% RS hit miss spike latency histogram
load(strcat(dataPath, 'RS.mat'))
sRS=find(strcmp({RS.subtype}, 'sRS'));

fs=16;
u2plot=sRS;
utitle='sRS';

yl=[0 0.0135];
figure('Position', [0 0 300 300])
hold all
for typi=1:numel(thrttalias)
  switch thrttalias{typi}
    case 'SubthreshHit'
      ttcol = [0.4 0.4 1];
      lw = .5;
    case 'SubthreshMiss'
      ttcol = [1 0.4 0.4];
      lw = .5;
    case 'SuprathreshHit'
      ttcol = [0 0 1];
      lw = 2.5;
    case 'SuprathreshMiss'
      ttcol = [1 0 0];
      lw = 2.5;
  end
  temp = [RS.FirstSpikeLatency];
  temp = [temp.(thrttalias{typi})];
  temp = [temp.PDF];
  
  plot(RS(1).FirstSpikeLatency.param.binnedSpikeLatencyNormProbtimeline, nanmean(temp(:,u2plot),2), 'Color', ttcol, 'LineWidth', lw)
end
xlim([0 50])
ylim(yl)
xlabel('1^s^t Spike Latency (ms)', 'FontSize', fs)
ylabel('Probability', 'FontSize', fs)
title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'))
set(gca, 'XTick', [0 50], 'YTick', 0:0.01:0.05, 'FontSize', fs)

text(0, yl(1)+.1*range(yl), ' Hit', 'FontSize', fs, 'Color', 'b', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom')
text(0, yl(1), ' Miss', 'FontSize', fs, 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom')
print('-painters', '-dpdf', strcat(figPath, 'fig8Ci'))


%%
u2plot=sRS;

tempReliability = struct();
for typi=1:numel(thrttalias)
temp = [RS.FirstSpikeLatency];
temp = [temp.(thrttalias{typi})];
tempReliability.(thrttalias{typi}) = [temp.Reliability];
end

yl=[0.005 0.032];
figure('Position', [0 0 150 300])
hold all
xt=[];
for tpi=1:2
  switch tpi
    case 1
      hpos=1;
      mpos=2;
      hp=tempReliability.SubthreshHit;
      mp=tempReliability.SubthreshMiss;
      hcol=[0.4 0.4 1];
      mcol=[1 0.4 0.4];
    case 2
      hpos=4;
      mpos=5;
      hp=tempReliability.SuprathreshHit;
      mp=tempReliability.SuprathreshMiss;
      hcol=[0 0 1];
      mcol=[1 0 0];
  end
  hcol='b';
  mcol='r';

  errorbar(hpos, nanmedian(hp(u2plot)), prctile(hp(u2plot),50)-prctile(hp(u2plot),25), prctile(hp(u2plot),75)-prctile(hp(u2plot),50), ...
    's-', 'Color', hcol, 'LineWidth', 2.5, 'MarkerFaceColor', hcol, 'MarkerSize',10, 'CapSize', 0)
  errorbar(mpos, nanmedian(mp(u2plot)), prctile(mp(u2plot),50)-prctile(mp(u2plot),25), prctile(mp(u2plot),75)-prctile(mp(u2plot),50), ...
    's-', 'Color', mcol, 'LineWidth', 2.5, 'MarkerFaceColor', mcol, 'MarkerSize',10, 'CapSize', 0)

  pval=signrank(hp(u2plot), mp(u2plot));
  if pval<0.05/2
    plot((hpos+mpos)/2, yl(2)-.09*range(yl), 'k+', 'LineWidth',1.5)
  end
  xt=[xt (hpos+mpos)/2];
end

if signrank(tempReliability.SubthreshHit(u2plot), tempReliability.SuprathreshHit(u2plot))<0.05/2
  plot([1 4], (yl(2)-.02*range(yl))*[1 1], 'k-', 'LineWidth',1)
  plot(2.5, yl(2), 'k*', 'LineWidth',1)
end
if signrank(tempReliability.SubthreshMiss(u2plot), tempReliability.SuprathreshMiss(u2plot))<0.05/2
  plot([2 5], (yl(2)-.06*range(yl))*[1 1], 'k-', 'LineWidth',1)
  plot(3.5, yl(2)-.04*range(yl), 'k*', 'LineWidth',1)
end

set(gca, 'YTick', 0:0.01:0.1, 'YTickLabelRotation', 90, 'XTick', xt, 'XTickLabel', {'Sub-', 'Supra-'}, 'XTickLabelRotation', 45, 'FontSize', fs, 'ActivePositionProperty', 'OuterPosition')
xlim([.5 5.5])
ylim(yl)
ylabel('1^s^t Spike Latency Reliability   ', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig8Cii'))

%% RS spike latency peak probability conditioned on FS CV
load(strcat(dataPath, 'RSFSpair.mat'))

temp = [RSFSpair.FS_Pre1sCVbins];
temp = [temp.Subthresh];
temp = [temp.RS_FirstSpikeLatency];
subfsrsvec = [temp.Reliability]';

temp = [RSFSpair.FS_Pre1sCVbins];
temp = [temp.Suprathresh];
temp = [temp.RS_FirstSpikeLatency];
supfsrsvec = [temp.Reliability]';

xl=[.5 3.5];
yl=[0.0045 0.013]; % 10ms bin
% yl=[0.0075 0.0165]; % 5ms bin

figure('Position', [0 0 900 300])
for spi=1:3
  switch spi
    case 1
      p2p=find(strcmp({RSFSpair.subtype_RS}, 'sRS') & strcmp({RSFSpair.subtype_FS}, 'sFS'));
      ptit='sFS-sRS pairs';
      fstitle='sFS';
    case 2
      p2p=find(strcmp({RSFSpair.subtype_RS}, 'sRS') & strcmp({RSFSpair.subtype_FS}, 'grnsFS'));
      ptit='grnsFS-sRS pairs';
      fstitle='grnsFS';
    case 3
      p2p=find(strcmp({RSFSpair.subtype_RS}, 'sRS') & strcmp({RSFSpair.subtype_FS}, 'plnsFS'));
      ptit='plnsFS-sRS pairs';
      fstitle='plnsFS';
  end
  
  subplot(1,3,spi)
  hold all
    
  plot(1:3, nanmean(subfsrsvec(p2p,2:4)), 'k-', 'LineWidth', 4)
  for ibin=2:4
%     submfc=(ibin-2)/3*[1 1 1];
    submfc=[0 0 0];
    errorbar(ibin-1, nanmean(subfsrsvec(p2p,ibin)), nanstd(subfsrsvec(p2p,ibin))/sqrt(numel(p2p)), ...
      'o-', 'Color', submfc, 'MarkerFaceColor', submfc, 'MarkerSize', 10, 'LineWidth', 4) %, 'CapSize', 0
  end
%   yl=ylim;
  
  tempmat=subfsrsvec(p2p, 2:end);
  [rnan, cnan]=find(isnan(tempmat));
  tempmat(unique(rnan),:)=[];
  
  [pfried, tblfried, statsfried]=friedman(tempmat, 1, 'off');
  c=multcompare(statsfried,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), yl(1)+.02*range(yl)+.08*range(yl)*[ci ci], 'k-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), yl(1)+.08*range(yl)*[ci ci], 'k*', 'LineWidth', 1)
  end
  text(xl(1)+0.05*range(xl), yl(1)+0.1*range(yl), sprintf('p=%.3f',pfried), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 14)
  text(xl(1)+0.05*range(xl), yl(1), sprintf('N=%d',numel(p2p)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold')
  
  
  plot(1:3, nanmean(supfsrsvec(p2p,2:4)), 'm-', 'LineWidth', 4)
  for ibin=2:4
    supmfc=[1 0 1];
    errorbar(ibin-1, nanmean(supfsrsvec(p2p,ibin)), nanstd(supfsrsvec(p2p,ibin))/sqrt(numel(p2p)), ...
      'o-', 'Color', supmfc, 'MarkerFaceColor', supmfc, 'MarkerSize', 10, 'LineWidth', 4) %, 'CapSize', 0
  end
  
  tempmat=supfsrsvec(p2p, 2:end);
  [rnan, cnan]=find(isnan(tempmat));
  tempmat(unique(rnan),:)=[];
  
  [pfried, tblfried, statsfried]=friedman(tempmat, 1, 'off');
  c=multcompare(statsfried,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), yl(2)+.02*range(yl)-.04*range(yl)*[ci ci], 'm-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), yl(2)+.04*range(yl)-.04*range(yl)*[ci ci], 'm*', 'LineWidth', 1)
  end
  text(xl(2), yl(2)-0.05*range(yl), sprintf('p=%.3f',pfried),'Color', 'm', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14)
  
  if spi==1
    text(xl(2), yl(1)+0.1*range(yl), 'Supra', 'Color', 'm', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
    text(xl(2), yl(1), 'Subthresh', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
  end
  
  xlim(xl)
  ylim(yl)
  set(gca, 'XTick', [1 3], 'XTickLabel', {'Low', 'High'}, 'FontSize', fs) % , 'YTick', [0.01 0.015]
  xlabel(strcat(fstitle,' CV Bin'), 'FontSize', fs)
  if spi==1
  ylabel('1^s^t Spike Latency Reliability', 'FontSize', fs)
  end
  %   title(strcat(ptit, ' N=', num2str(numel(p2p))))
  
end
print('-painters', '-dpdf', strcat(figPath, 'fig8E'))

%%
gausssigma=3;
gausswindow=21;
tempfil= exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter=tempfil/sum(tempfil);

hlab=' Hit';
hcol='b';
mlab=' Miss';
mcol='r';

temp = [RS.PSTH];
suphsp=conv2([temp.SuprathreshHit], tempfilter, 'same');
supmsp=conv2([temp.SuprathreshMiss], tempfilter, 'same');

subhsp=conv2([temp.SubthreshHit], tempfilter, 'same');
submsp=conv2([temp.SubthreshMiss], tempfilter, 'same');

u2plot=sRS;
utitle='sRS';
yl=[4 18];
figure('Position', [100 100 300 300])
hold all
plot(RS(1).PSTH.timeline, nanmean(suphsp(:,u2plot),2),'Color', [0 0 1], 'LineWidth', 2.5)
plot(RS(1).PSTH.timeline, nanmean(supmsp(:,u2plot),2),'Color', [1 0 0], 'LineWidth', 2.5)
plot(RS(1).PSTH.timeline, nanmean(subhsp(:,u2plot),2),'Color', [.4 .4 1], 'LineWidth', .5)
plot(RS(1).PSTH.timeline, nanmean(submsp(:,u2plot),2),'Color', [1 .4 .4], 'LineWidth', .5)
xlim([-250 100])
ylim(yl)
% title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)
set(gca, 'XTick', [-250 -125 0 100 500], 'FontSize', fs)
ylabel('Rate (Hz)', 'FontSize', fs)
text(-250,yl(2), hlab, 'FontSize', 16, 'FontWeight', 'bold', 'Color', hcol, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(-250,yl(2)-.1*range(yl), mlab, 'FontSize', 16, 'FontWeight', 'bold', 'Color', mcol, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
xlabel('Time from Sensory Onset (ms)     ', 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig8Di'))


%%
u2plot=sRS;

tsttoi= RS(1).PSTH.timeline>=0 & RS(1).PSTH.timeline<=100;
tsttoiname= 'Post 100ms';
yl=[0 25];

figure('Position', [0 0 150 300])
hold all
xt=[];
for tpi=1:2
  switch tpi
    case 1
      hpos=1;
      mpos=2;
      hp=mean(subhsp(tsttoi,:),1);
      mp=mean(submsp(tsttoi,:),1);
    case 2
      hpos=4;
      mpos=5;
      hp=mean(suphsp(tsttoi,:),1);
      mp=mean(supmsp(tsttoi,:),1);
  end
  hcol='b';
  mcol='r';
  
  errorbar(hpos, nanmedian(hp(u2plot)), prctile(hp(u2plot),50)-prctile(hp(u2plot),25), prctile(hp(u2plot),75)-prctile(hp(u2plot),50), ...
    's-', 'Color', hcol, 'LineWidth', 2.5, 'MarkerFaceColor', hcol, 'MarkerSize',10, 'CapSize', 0)
  errorbar(mpos, nanmedian(mp(u2plot)), prctile(mp(u2plot),50)-prctile(mp(u2plot),25), prctile(mp(u2plot),75)-prctile(mp(u2plot),50), ...
    's-', 'Color', mcol, 'LineWidth', 2.5, 'MarkerFaceColor', mcol, 'MarkerSize',10, 'CapSize', 0)
  
  pval=signrank(hp(u2plot), mp(u2plot));
  if pval<0.05/2
  plot((hpos+mpos)/2, yl(2)-.1*range(yl), 'k+', 'LineWidth',1.5)
  end
  xt=[xt (hpos+mpos)/2];
end

if signrank(mean(subhsp(tsttoi,u2plot),1), mean(suphsp(tsttoi,u2plot),1))<0.05/2
  plot([1 4], (yl(2)-.02*range(yl))*[1 1], 'k-', 'LineWidth',1)
  plot(2.5, yl(2), 'k*', 'LineWidth',1)
end
if signrank(mean(submsp(tsttoi,u2plot),1), mean(supmsp(tsttoi,u2plot),1))<0.05/2
  plot([2 5], (yl(2)-.06*range(yl))*[1 1], 'k-', 'LineWidth',1)
  plot(3.5, yl(2)-.04*range(yl), 'k*', 'LineWidth',1)
end

set(gca, 'YTickLabelRotation', 90, 'XTick', xt, 'XTickLabel', {'Sub-', 'Supra-'}, 'XTickLabelRotation', 45, 'FontSize', fs, 'ActivePositionProperty', 'OuterPosition')
xlim([.5 5.5])
ylim(yl)
ylabel(strcat(tsttoiname, ' Rate (Hz)'), 'FontSize', fs)
print('-painters', '-dpdf', strcat(figPath, 'fig8Dii'))


%% CVsq bins spike count
tsttoi= RSFSpair(1).FS_Pre1sCVbins.param.PSTHTimeline>0 & RSFSpair(1).FS_Pre1sCVbins.param.PSTHTimeline<=100;
tsttoiname= 'Post 100ms';

temp = [RSFSpair.FS_Pre1sCVbins];
temp = [temp.Subthresh];
submat = cat(3, temp.RS_PSTH);

temp = [RSFSpair.FS_Pre1sCVbins];
temp = [temp.Suprathresh];
supmat = cat(3, temp.RS_PSTH);

subfsrsvec = squeeze(mean(submat(tsttoi,:,:),1))';
supfsrsvec = squeeze(mean(supmat(tsttoi,:,:),1))';

figure('Position', [0 0 900 300])
for spi=1:3
  switch spi
    case 1
      p2p=find(strcmp({RSFSpair.subtype_RS}, 'sRS') & strcmp({RSFSpair.subtype_FS}, 'sFS'));
      ptit='sFS-sRS pairs';
      fstitle='sFS';
      yl=[2 11];
    case 2
      p2p=find(strcmp({RSFSpair.subtype_RS}, 'sRS') & strcmp({RSFSpair.subtype_FS}, 'grnsFS'));
      ptit='grnsFS-sRS pairs';
      fstitle='grnsFS';
      yl=[1 6];
    case 3
      p2p=find(strcmp({RSFSpair.subtype_RS}, 'sRS') & strcmp({RSFSpair.subtype_FS}, 'plnsFS'));
      ptit='plnsFS-sRS pairs';
      fstitle='plnsFS';
      yl=[2 11];
  end
  
  subplot(1,3,spi)
  hold all
    
  plot(1:3, nanmean(subfsrsvec(p2p,2:4)), 'k-', 'LineWidth', 4)
  for ibin=2:4
%     submfc=(ibin-2)/3*[1 1 1];
    submfc=[0 0 0];
    errorbar(ibin-1, nanmean(subfsrsvec(p2p,ibin)), nanstd(subfsrsvec(p2p,ibin))/sqrt(numel(p2p)), ...
      'o-', 'Color', submfc, 'MarkerFaceColor', submfc, 'MarkerSize', 10, 'LineWidth', 4) %, 'CapSize', 0
  end
%   yl=ylim;
  
  tempmat=subfsrsvec(p2p, 2:end);
  [rnan, cnan]=find(isnan(tempmat));
  tempmat(unique(rnan),:)=[];
  
  [pfried, tblfried, statsfried]=friedman(tempmat, 1, 'off');
  % c=multcompare(statsfried,'CType','dunn-sidak','Display','off');
  c=multcompare(statsfried,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), yl(1)+.02*range(yl)+.08*range(yl)*[ci ci], 'k-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), yl(1)+.08*range(yl)*[ci ci], 'k*', 'LineWidth', 1)
  end
  text(xl(1)+0.05*range(xl), yl(1)+0.1*range(yl), sprintf('p=%.3f',pfried), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 14)
  text(xl(1)+0.05*range(xl), yl(1), sprintf('N=%d',numel(p2p)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold')

  
  plot(1:3, nanmean(supfsrsvec(p2p,2:4)), 'm-', 'LineWidth', 4)
  for ibin=2:4
    supmfc=[1 0 1];
    errorbar(ibin-1, nanmean(supfsrsvec(p2p,ibin)), nanstd(supfsrsvec(p2p,ibin))/sqrt(numel(p2p)), ...
      'o-', 'Color', supmfc, 'MarkerFaceColor', supmfc, 'MarkerSize', 10, 'LineWidth', 4) %, 'CapSize', 0
  end
  
  tempmat=supfsrsvec(p2p, 2:end);
  [rnan, cnan]=find(isnan(tempmat));
  tempmat(unique(rnan),:)=[];
  
  [pfried, tblfried, statsfried]=friedman(tempmat, 1, 'off');
  c=multcompare(statsfried,'CType','bonferroni','Display','off');
  sigp=find(c(:,6)<0.05);
  for ci=1:numel(sigp)
    plot(c(sigp(ci),1:2), yl(2)+.02*range(yl)-.04*range(yl)*[ci ci], 'm-', 'LineWidth', 1)
    plot(mean(c(sigp(ci),1:2)), yl(2)+.04*range(yl)-.04*range(yl)*[ci ci], 'm*', 'LineWidth', 1)
  end
  text(xl(2), yl(2)-0.05*range(yl), sprintf('p=%.3f',pfried),'Color', 'm', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14)

  if spi==1
    text(xl(2), yl(1)+0.1*range(yl), 'Supra', 'Color', 'm', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
    text(xl(2), yl(1), 'Subthresh', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
  end
  
  xl=[.5 3.5];
  xlim(xl)
  ylim(yl)
  set(gca, 'XTick', [1 3], 'XTickLabel', {'Low', 'High'}, 'FontSize', fs) % , 'YTick', [0.01 0.015]
  xlabel(strcat(fstitle,' CV Bin'), 'FontSize', fs)
  if spi==1
  ylabel(strcat(tsttoiname, ' Rate (Hz)'), 'FontSize', fs)
  end
  %   title(strcat(ptit, ' N=', num2str(numel(p2p))))  
end
print('-painters', '-dpdf', strcat(figPath, 'fig8F'))
