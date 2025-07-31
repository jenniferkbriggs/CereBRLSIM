%% This script uses data from - which can be obtain via reasonable request to Philip Ainslie: 
% Carr, Jay MJR, et al. "Internal carotid and brachial artery 
% shear‐dependent vasodilator function in young healthy humans." 
% The Journal of physiology 598.23 (2020): 5333-5350. Figure 2B

j = 1; %for indexing if you are going through more than 1 data sets
endot= [20]; %endothelial response time (tau_endo)
cuttime = [150]
dat = readtable(['cFMD-Table 1.csv'], 'ReadRowNames',false, 'ReadVariableNames', true);

CBFv1 = dat.ShearRate_1_s_;
diameter = dat.ContinuousDiameter_cm_;
time = dat.Time;


%Smooth and interpolate: 
CBFv = movmean(CBFv1, 10, 'omitmissing');
diam = movmean(diameter, 10, 'omitmissing');
CBF(:,1) = time;

cPerc = (CBFv)/CBFv(1)-1;
dPerc = diam/diam(1)-1;

%convert shear stress to CBF
%convert diam to radius of 0.1: 
sc = diam(1,1)/0.1;
CBF(:,2) = CBFv.*pi.*(diam./sc).^2;%Divide by 4 to convert to radius of 0.1
CBF(:,2) = CBF(:,2);%+0.6*rand(size(CBF(:,2)));

%remove nans: 
cut = find(isnan(CBF(:,2)));
time(cut) = [];
CBF(cut,:) = [];
diam(cut) = [];
CBFv(cut) = [];

cut = find(isnan(CBF(:,1)));
time(cut) = [];
CBF(cut,:) = [];
diam(cut) = [];
CBFv(cut) = [];

cut = find(time > cuttime(j))
time(cut) = [];
CBF(cut,:) = [];
diam(cut) = [];
CBFv(cut) = [];

cut = find(time < cutstart(j))
time(cut) = [];
CBF(cut,:) = [];
diam(cut) = [];
CBFv(cut) = [];


Press = 70.*ones(size(CBF)); %Assume pressure is normal CPP = 70 mmHg

load("FakeData/FMD/Paramvals.mat")

%% 
paramvals(7) = 1; %endothelial mechanism
paramvals(4) = endot(j);
paramvals(8) = 1;
paramvals(6) = 1;
paramvals(5) = 1;
paramvals(2) = CBF(1,2);
try
    paramvals(10) = ECO2(1,2);
catch
    paramvals(10) = 40;
end
% 
% paramvals(7) = 1;
% paramvals(1) = 0.1;
% paramvals(3) =1.5;
% paramvals(14) = 6;
% paramvals(4) = 35


IC = [paramvals(2), paramvals(1), 0, 0, 0, 0, 0, 0, 0];
%figure(2), plot(time, CBF(:,2)), ylabel('Flow'), yyaxis right, plot(time, diam), ylabel('Diam')
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
        paramvals(6:8) = [0,1,1]
        tit ='Endothelial and Metabolic Mechanisms'
        ct = 3
    elseif i == 4
        paramvals(6:8) = [1,1,0]
        tit ='Endothelial and Myogenic Mechanisms'
        ct = 4
    elseif i == 5
        paramvals(6:8) = [0,1,0]
        tit ='Endothelial Mechanism Only'
        ct = 2
    end
[t,y] = ode23(@(t, y) cerebralbloodflow_changeshearstress(t,y,paramvals,[time, Press/75], CBF, [], [], []), time(4:end)', IC);

%scale diameters:
 if i == 1
b1 = [ones(size(y(:,2))) y(:,2)]\diam(4:end);

%or do just standard scaling: 
%b1 = [1;1];

if length(find(isnan(b1))>1)
    b1 = [ones(size(y(:,2))) y(:,2)]\diam(1:end-3);
end
 end
sim_diam = [ones(size(y(:,2))) y(:,2)]*(b1);%if it's really bad sometimes the slope is -

if 1 %turn on to plot
fig = figure, nexttile,
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
ylabel('Diameter (mm)')

set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
xlabel('Time (s)')

title(tit)
if i == 1
    title(["Metabolic, Myogenic,"; "and Endothelial Mechanisms"])
    ax = gca
    ylims = ax.YLim
elseif i == 2
    legend('Scaled Diameter_{sim}','Diameter_{data}');%, 'CBF_{data}')%, title(['Radius = ', num2str(paramvals(1))])
end
ylim(ylims);


if j == 1
    if i == 3 || i == 5
        axes('Position',[0.55, 0.6, 0.35, 0.3])
        plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
        ylabel('Diameter')
        xlabel('              Time')
        set(gca, 'XTick', [])
        set(gca, 'YColor', 'k')
        set(gca, 'box','off')
    end
elseif j == 2
    if i == 3 || i == 5
        axes('Position', [0.6500 0.2000 0.2500 0.3000])
        plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
        ylabel('Diameter')
        xlabel('Time')
        set(gca, 'XTick', [])
        set(gca, 'YColor', 'k')
        set(gca, 'box','off')
    end
elseif j == 3
    if i == 3 || i == 5
        axes('Position', [0.2200 0.6400 0.2500 0.2500])
        yyaxis right
        plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
        ylabel('Diameter')
        xlabel('Time')
        set(gca, 'XTick', [])
        set(gca, 'YColor', 'k')
        set(gca, 'box','off')
        yyaxis left
        set(gca, 'YColor', 'none')
    end
elseif j == 4
    if i == 3 || i == 5
        axes('Position', [0.6500 0.2000 0.2500 0.2500])
        plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
        ylabel('Diameter')
        xlabel('Time')
        set(gca, 'XTick', [])
        set(gca, 'YColor', 'k')
        set(gca, 'box','off')
    end
elseif j == 5
    if i == 3 || i == 5
        axes('Position', [0.2500 0.6400 0.2500 0.2500])
        yyaxis right
        plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
        ylabel('Diameter')
        xlabel('Time')
        set(gca, 'XTick', [])
        set(gca, 'YColor', 'k')
        set(gca, 'box','off')
        yyaxis left
        set(gca, 'YColor', 'none')
    end
end
% saveas(gcf, [filename, '_', strrep(tit, ' ', ''), '.fig'])
% saveas(gcf, [filename, '_', strrep(tit, ' ', ''), '.png'])


%plot Radius
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
yyaxis right, plot(t-t(1), y(:,3:5), 'linewidth',3)
if i == 2
l = legend('\phi_{myo}', '\phi_{endo}', '\phi_{meta}')%, title(['Radius = ', num2str(paramvals(1))])
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
end
yyaxis left
ax2 = gca
set(ax2.YAxis(1), 'Visible', 'off')
yyaxis right
ylim(yylimR)
ylabel('Force From Mechansims')
% saveas(gcf, [filename, '_', strrep(tit, ' ', ''), 'Radius.fig'])
% saveas(gcf, [filename, '_', strrep(tit, ' ', ''), 'Radius.png'])


if i == 1
%plot pressure: 
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot(time, CBFv, 'Color', [190, 94, 121]./215, 'linewidth',3)
ylabel('Shear Rate (1/s)')
xlabel('Time (s)')
set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
% saveas(fig, [filename, '_Shear.fig'])
% saveas(fig, [filename, '_Shear.png'])
end
end

corelation = corr(sim_diam, diam(4:end))
mse = mean((sim_diam-diam(4:end)).^2)


mse_all(ct) = [mse];
corr_all(ct)= [corelation];

end
