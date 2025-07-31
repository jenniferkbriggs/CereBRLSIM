%Run Experiment to replicate Mahony and Evans:
clear all
close all
clc
try
    cd('/Users/jkbriggs/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Anschutz/Albers/CVTR_Model/TympkoData/CSV_data/NVC')
catch
    cd('~/OneDrive - The University of Colorado Denver/Anschutz/Albers/CVTR_Model/TympkoData/CSV_data/NVC')
end



filename = '02'
Flow = readmatrix(['FakeData/NVC/MCAv.csv']); %Row 1 is time, Row 2 is data
markers = readmatrix(['FakeData/NVC/timestamps.csv']); %time when experiment starts, eyes open, and eyes close 
ABP = readmatrix(['FakeData/NVC/ABP.csv']);  %Row 1 is time, Row 2 is data


%interpolate data
t = ABP(:,1)';
y = ABP(:,2)'; %assume steady icp

%cut data to first marker (e.g. time when experiment starts, eyes open, and eyes close): 
y(t<markers(1,1)) = [];
t(t<markers(1,1)) = [];

Flow(Flow(:,1)<markers(1,1),:) = [];


%clean up data by sorting: 
[t2, i] = sort(t);
y = y(i);
n=unique(t2);
if length(n) < length(t2)
    indx = find(diff(t2)==0)
    t2(indx) = [];
    y(indx) = [];
end
t = t2;


Flow_old = Flow;
%clean up data by sorting: 
tflow = Flow(:,1);
[t2, i] = sort(tflow);
Flow2 = Flow(i,2);
n=unique(t2);
if length(n) < length(t2)
    indx = find(diff(t2)==0)
    t2(indx) = [];
    Flow2(indx) = [];
end
Flow = [t2, Flow2];


%interpolate flow:
Datasmoothed = movmean(Flow(:,2), 5);
Flow = interp1(Flow(:,1),Flow(:,2), t)% + rand(size(t));

%remove nans
nanindx = find(isnan(Flow));
t(nanindx) = [];
y(nanindx) = [];
Flow(nanindx) = [];

%% -- calculate metabolic change -- % 
%find time that eyes are closed and calculate lagged cross correlation:
oo = markers(2,1);
oo_p = find(t<oo);
oo_p = oo_p(end);
res = mean(diff(t))
[x,lag] = xcov(Flow(1:oo_p), y(1:oo_p), round(80/res)); %max lag is 10 seconds
[mx, indx] = max(x);
lag_s = lag(indx)*res; %this is the lag in seconds: 

%Metabolic change - increase when eyes open and decrease when eyes close
P = @(t, tau,mi,ma,theta) rescale((1./(1+exp(theta*(t-tau)))),mi,ma);
theta = 2; 
tau = markers(2,1)+2; %start of the bump
mi = 0; %min value of P in (mmHg)
ma = 1.5; %max value = ma + mi (mmHg)
t2 = [t(1):mean(diff(t)):t(end)]; %make this smoother
qn_nodiff = (P(t2, tau,mi,ma+mi, -theta)-P(t2, markers(3,1)+6,mi,ma+mi, -theta/5));%
M = diff(qn_nodiff).*250;%

IC = [paramvals(2), paramvals(1), 0, 0, 0, 0, 0, 0, 0];
paramvals(5) = 1;

paramvals(6) = 0;
paramvals(7) = 0;
paramvals(8) = 1;

for i = 1:5
    if i == 1
        paramvals(6:8) = [1,1,1]
        tit ='Metabolic, Myogenic, and Endothelial Mechanisms'
        ct = 5
    elseif i == 2
        paramvals(6:8) = [0,0,0]
        tit ='No CVTR'
        ct = 1
    elseif i == 3
        paramvals(6:8) = [1,0,1]
        tit ='Metabolic and Myogenic Mechanisms'
        ct = 3
    elseif i == 4
        paramvals(6:8) = [0,1,1]
        tit ='Metabolic and Endothelial Mechanisms'
        ct = 4
    elseif i == 5
        paramvals(6:8) = [0,0,1]
        tit ='Metabolic Mechanism Only'
        ct = 2
    end

%NOTE: requires function to change blood flow demand in the o
[time,state] = ode23(@(time, state) CereBRLSIM_NVC(time, state,paramvals,[t', y'/75], [t2(2:end)' qn']), t(4:end)', IC);

%find the scaling between cbfv and simulated flow:
x = 0.1
rad_perc = (squeeze(state(:,2)))./squeeze(state(1,2)); %percent change of the radius
rad_perc = rad_perc-1; %quantify if it is getting bigger or smaller.  - this is what we will rescale for the MCA 
MCA_rad = (rad_perc.*x+1);
sim_cbfv = squeeze(state(:,1,:))./(MCA_rad.^2*pi);



%if i == 1 %scale based on metabolic only
    mdl = fitlm(smoothdata(sim_cbfv, 'gaussian'), smoothdata(Flow(4:end), 'gaussian'));
%end

if 0
%plot model vs data
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot((time-time(1)), mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*sim_cbfv, 'k:', 'linewidth',3), hold on,
plot((t-t(1)), Flow(1:end), 'k', 'linewidth',3)
if i == 2
legend('Scaled CBv_{sim}', 'CBv_{data}')%, title(['Radius = ', num2str(paramvals(1))])
end
set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
xlabel('Time (s)')
ylabel('CBv (cm/s)')
%tit = input("Input Title and Save Name   ", "s")
title(tit)
if i == 1
    title(["Metabolic, Myogenic,"; "and Endothelial Mechanisms"])
end
saveas(gcf, [filename, '_', strrep(tit, ' ', ''), '.fig'])
saveas(gcf, [filename, '_', strrep(tit, ' ', ''), '.png'])



%plot model vs data
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot((time-time(1)), state(:,2), 'k', 'linewidth',3)
yyaxis right, plot(time-time(1), state(:,3:5), 'linewidth',3)
if i == 2
l = legend('Radius', '\phi_{myo}', '\phi_{endo}', '\phi_{meta}')%, title(['Radius = ', num2str(paramvals(1))])
%keyboard
end
set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
xlabel('Time (s)')
%tit = input("Input Title and Save Name   ", "s")
title(tit)
if i == 1
    title(["Metabolic, Myogenic,"; "and Endothelial Mechanisms"])
    yyaxis right
    ax = gca
    yylimR = ax.YLim
    yyaxis left
    ax = gca
    yylimL = ax.YLim
end
yyaxis left
ylim(yylimL)
ylabel('Simulated Radius')
yyaxis right
ylim(yylimR)
ylabel('Force From Mechansims')
saveas(gcf, [filename, '_', strrep(tit, ' ', ''), 'Radius.fig'])
saveas(gcf, [filename, '_', strrep(tit, ' ', ''), 'Radius.png'])


if i == 1
%plot pressure: 
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot(t-t(1), y, 'Color', [190, 94, 121]./215, 'linewidth',3)
ylabel('ABP (mmHg)')
xlabel('Time (s)')
set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
saveas(gcf, [filename, '_ABP.fig'])
saveas(gcf, [filename, '_ABP.png'])

%plot Metabolism: 
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
%plot(t(2:end)-t(1), qn, 'Color', [190, 94, 121]./215, 'linewidth',3, 'LineStyle',":")
plot(t2-t2(1), qn_nodiff, 'Color', [190, 94, 121]./215, 'linewidth',3)
xline(markers(2,1)-t(1), 'label', 'Eyes Open', 'FontSize',15, 'LineWidth',2, 'LabelVerticalAlignment','middle', 'Color','k','LineStyle',":", 'LabelHorizontalAlignment','left')
xline(markers(3,1)-t(1), 'label', 'Eyes Closed', 'FontSize',15, 'LineWidth',2, 'LabelVerticalAlignment','middle', 'Color','k','LineStyle',":",'LabelHorizontalAlignment','left')
ylabel('$\mathcal{M}$', 'Interpreter','latex')
set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
saveas(gcf, [filename, '_M.fig'])
saveas(gcf, [filename, '_M.png'])

end
end


corelation = corr(mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*sim_cbfv, Flow(4:end)')
mse = mean((mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*sim_cbfv-Flow(4:end)').^2);
avg_pressure = mean(y)
range_pressure = range(y)
avg_cbv = mean(Flow)
range_cbv = range(Flow)

end
