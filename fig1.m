figPath = uigetdir('', 'Select folder for figures to be saved into');
figPath = strcat(figPath, '\');
disp(figPath)

dataPath = uigetdir('', 'Select folder where data was downloaded to');
dataPath = strcat(dataPath, '\');
disp(dataPath)

load(strcat(dataPath, 'FS.mat'))

sFS=find(strcmp({FS.subtype}, 'sFS'));
plnsFS=find(strcmp({FS.subtype}, 'plnsFS'));
grnsFS=find(strcmp({FS.subtype}, 'grnsFS'));
nsFS=find(contains({FS.subtype}, 'nsFS'));
supFS=find(strcmp({FS.subtype}, 'suppressedFS'));

%% trial structure schematic
fs=14;
load(strcat(dataPath, 'HSVGATDR_151025_trial_52.mat'));
tt=[1:size(data,1)]/26000 -4.15;
rw=ones(1,size(data,1));
rw(tt==0)=0;
rw(tt>0 & tt<=.7)=0.6;
yh=3.9;
figure('Position', [300 300 400 120])
hold all
% imagesc(tt, [.9 2.7], repmat(rw,2,1),'AlphaData',0.8)
% colormap hot
plot([0 .7], [yh yh], '-', 'Color',[1 0.7 0], 'LineWidth', 3)
plot(tt, -data(:,2), 'k-', 'LineWidth', 1)
text(.35,yh+.15,'Report Window', 'Color',[1 0.7 0], 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
plot([.9 1], [-.8 -.8], 'k-', 'LineWidth', 2)
text(.95,-.8,'100ms', 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
plot([0 0], [-1 -.4], 'r-', 'LineWidth', 2)
plot(0, -.4, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
text(0,-.8,'Sensory Onset', 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
plot([.75 .75], [0 3.6], 'm-', 'LineWidth', 2)
plot(.75,3.2, 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 6)
plot(.75,0.4, 'mv', 'MarkerFaceColor', 'm', 'MarkerSize', 6)
text(.8,1.8,{'Stimulus', 'Amplitude'}, 'Color','m', 'FontSize', fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
text(-.540,1.8,{'Vibrissae', 'Stimulus'}, 'FontWeight', 'bold', 'Color','k', 'FontSize', fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
xlim([-.525 1.2])
ylim([-1 yh])
axis off
title(' ', 'FontSize', fs)

print('-painters', '-dpdf', strcat(figPath, 'fig1A'))

%%
fs=14; pltopt=0;
gausssigma=3;
gausswindow=21;
tempfil= exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter=tempfil/sum(tempfil);

TSTtimeline = FS(1).PSTH.timeline;
toi=TSTtimeline>-475 & TSTtimeline<=1000;

temp = [FS.PSTH];

figure('Position', [300 300 400 120])
hold all
for sa=1:4
  switch sa
    case 1
      tempst=[temp.Catch];
    case 2
      tempst=[temp.Sub];
    case 3
      tempst=[temp.Supra];
    case 4
      tempst=[temp.Maxstim];
  end
  tempst=conv(mean(tempst(:,sFS),2), tempfilter, 'same');
  st2p=tempst(toi);
  plot(TSTtimeline(toi), st2p, 'Color', (sa-1)/3*[1 0 1], 'LineWidth', 1);
  plot([1050 1050], sa*6+[2 8], '-', 'Color', (sa-1)/3*[1 0 1], 'LineWidth', 5)
end
xlim(1000*[-.525 1.2])
axis off
ylim([8 32])
plot([0 0], [8 10], 'r-', 'LineWidth',2)
text(-540,20,{'sFS', 'Spike Rate'}, 'FontWeight', 'bold', 'Color','k', 'FontSize', fs, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
h1=text(1075, 20, 'Stim Amp', 'FontSize', fs, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
set(h1, 'Rotation', -90)
text(0, 9, '0', 'FontSize', fs, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')
if pltopt==1
plot([-525 -500], [10 10], 'k-', 'LineWidth',2)
plot([-525 -500], [30 30], 'k-', 'LineWidth',2)
text(-525, 10, '10', 'FontSize', fs, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
text(-525, 30, '30', 'FontSize', fs, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
h2=text(-575,20, 'Hz', 'FontSize', fs, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
set(h2, 'Rotation', 90)
else
% plot([750 750], [15 25], 'k-', 'LineWidth', 2)
errorbar(750, 20, 5, 'k-', 'LineWidth', 2)
text(800, 15.5, '15', 'FontSize', fs, 'VerticalAlignment', 'middle')
text(800, 25.5, '25', 'FontSize', fs, 'VerticalAlignment', 'middle')
text(800, 20.5, 'Hz', 'FontSize', fs, 'FontAngle', 'italic', 'FontWeight', 'bold', 'VerticalAlignment', 'middle')
end

print('-painters', '-dpdf', strcat(figPath, 'fig1Aa'))

%% SP histogram
fs=16;

temp=[FS.SP];
temp=[temp.Post100ms];
temp=[temp.Rate];

yl=[0 60];
figure('Position', [0 0 300 300])
hold all
for spi=1:3
  switch spi
    case 3
      u2plot=supFS;
      uc=[0 0 0];
      utitle='Suppressed FS';
    case 2
      u2plot=nsFS;
      uc=[0 .5 0];
      utitle='Non-Sensory FS';
    case 1
      u2plot=sFS;
      uc=[.5 .5 .5];
      utitle='Sensory FS';
  end

  if spi==1
    [hc,he]=histcounts(temp(u2plot), 0:0.04:1);
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(temp(u2plot), he, 'FaceColor', uc, 'FaceAlpha', .5, 'EdgeColor', 'none')
  end
text(he(2),yl(2)-.1*range(yl)*(spi-1), utitle, 'FontSize', 16, 'FontWeight', 'bold', 'Color', uc(1:3), 'VerticalAlignment', 'bottom')
end
set(gca, 'YTick', 0:20:40, 'FontSize', fs)
ylim(yl)
xlabel('Stimulus Probability', 'FontSize', fs)
ylabel('Count', 'FontSize', fs)

print('-painters', '-dpdf', strcat(figPath, 'fig1B'))

%% neurometric spike count FS
temp=[FS.SpikeRate];
temp=[temp.Post100ms];
thrmat=[[temp.Catch]; [temp.Sub]; [temp.Supra]; [temp.Maxstim]];

lw=4;
yl=[10 25];

for pltopt=1:2

switch pltopt
  case 1
    u2plot=sFS;
    utitle='sFS';
  case 2
    u2plot=nsFS;
    utitle='nsFS';
  case 3
    u2plot=supFS;
    utitle='Suppressed FS';
end

figure('Position', [0 0 270 300])
hold all
plot(1:size(thrmat,1), nanmean(thrmat(:,u2plot),2), 'k-', 'LineWidth', lw)
for ii=1:size(thrmat,1)
errorbar(ii, nanmean(thrmat(ii,u2plot),2), nanstd(thrmat(ii,u2plot),0,2)/sqrt(numel(u2plot)), 'o-', 'Color', (ii-1)/(size(thrmat,1)-1)*[1 0 1], 'MarkerFaceColor', (ii-1)/(size(thrmat,1)-1)*[1 0 1], 'LineWidth', lw, 'MarkerSize', 8)
end
ylim(yl)
xlim([0.5 size(thrmat,1)+0.5])
set(gca, 'XTick', [1 size(thrmat,1)], 'XTickLabel', {'Catch', 'Max'}, 'FontSize', fs)
xlabel('Stimulus Amplitude')
ylabel('Rate (Hz)', 'FontSize', fs)
title(strcat(utitle, ' (N=', num2str(numel(u2plot)), ')'), 'FontSize', fs)

print('-painters', '-dpdf', strcat(figPath, 'fig1C_', utitle))
end

%% spike latency histogram
temp=[FS.FirstSpikeLatency];
temp=[temp.Maxstim];
temp=[temp.Timing];

he=0:3:250;
xl=[3 40];
figure('Position', [0 0 300 300])
hold all
for spi=1:3
  switch spi
    case 3
      u2plot=supFS;
      uc=[0 0 0];
      utitle='supFS';
    case 2
      u2plot=nsFS;
      uc=[0 .5 0];
      utitle='nsFS';
    case 1
      u2plot=sFS;
      uc=[.5 .5 .5];
      utitle='sFS';
  end
  if spi==1
    [hc,he]=histcounts(temp(u2plot), 0:3:250); 
    stairs(he(1:end-1), hc, 'Color', uc, 'LineWidth',4)
  else
    histogram(temp(u2plot), he, 'FaceColor', uc, 'FaceAlpha', .5, 'EdgeColor', 'none')
  end
end
set(gca, 'XTick', 0:10:100, 'YTick', 0:10:40, 'FontSize', fs)
xlim(xl)
xlabel('1^s^t Spike Latency (ms)', 'FontSize', fs)
ylabel('Count', 'FontSize', fs)

print('-painters', '-dpdf', strcat(figPath, 'fig1D'))

%% spike latency distributions
temp=[FS.FirstSpikeLatency];
thrstimtypes={'Catch', 'Sub', 'Supra', 'Maxstim'};
thrmat=[];
for typi=1:numel(thrstimtypes)
  temp1 = [temp.(thrstimtypes{typi})];
  thrmat = cat(3, thrmat, [temp1.PDF]);
end

for pltopt=1:2

switch pltopt
  case 1
    u2plot=sFS;
    utitle='sFS';
    bv='v';
  case 2
    u2plot=nsFS;
    utitle='nsFS';
  case 3
    u2plot=supFS;
    utitle='supFS';
end


figure('Position', [pltopt*300 100 300 300])
hold all
for ii=1:size(thrmat,3)
  plot(FS(1).FirstSpikeLatency.param.binnedSpikeLatencyNormProbtimeline, mean(thrmat(:,u2plot,ii),2), 'Color', (ii-1)/(size(thrmat,3)-1)*[1 0 1], 'LineWidth', 2.5)
  if pltopt==1
    switch bv
      case 'h' %horizontal
        plot(5*ii+[25 30], [0.025 0.025], '-', 'Color', (ii-1)/(size(thrmat,3)-1)*[1 0 1], 'LineWidth', 5)
        if ii==size(thrmat,3)
          plot(49, 0.025, '>', 'MarkerSize', 10, 'MarkerFaceColor', (ii-1)/(size(thrmat,3)-1)*[1 0 1], 'Color', (ii-1)/(size(thrmat,3)-1)*[1 0 1])
        end
      case 'v' %vertical
        plot([43 43], 0.002*ii+[0.015 0.017], '-', 'Color', (ii-1)/(size(thrmat,3)-1)*[1 0 1], 'LineWidth', 5)
    end
  end
end
if pltopt==1
  switch bv
    case 'h' %horizontal
      text(30, 0.025, 'Stim Amp', 'FontSize', 14, 'VerticalAlignment', 'top')
    case 'v' %vertical
      h1=text(50, 0.021, {'Stim Amp'}, 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
      set(h1, 'Rotation', -90)
  end
end
xlim([0 50])
yl=[0 0.025];
ylim(yl)
xlabel('1^s^t Spike Latency (ms)', 'FontSize', fs);
ylabel('Probability', 'FontSize', fs)
set(gca, 'XTick', 0:25:100, 'YTick', yl,'FontSize', fs)%, 'ActivePositionProperty', 'outerposition')
ax = gca();
ax.YRuler.Exponent = -2;

print('-painters', '-dpdf', strcat(figPath, 'fig1E_', utitle))
end
