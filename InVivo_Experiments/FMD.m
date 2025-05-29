%% This script uses data extracted from: 
% Carr, Jay MJR, et al. "Internal carotid and brachial artery 
% shear‐dependent vasodilator function in young healthy humans." 
% The Journal of physiology 598.23 (2020): 5333-5350. Figure 2B

% Since we are driving with shear stress rather than pressure, we will use
% the model itteration called "cerebralbloodflow_changeshearstress.m".

% ------- load data ------- 
CBFv1 = readmatrix('CBFv_Hoiland.csv');
diameter = readmatrix('ICAdiam_Hoiland.csv');
time = [CBFv1(1,1):0.5:CBFv1(end,1)]';

%Interpolate true values. 
CBFv = interp1(CBFv1(:,1), CBFv1(:,2), time, 'linear');
diam = interp1(diameter(:,1), diameter(:,2), time, 'linear');

%calculate shear stress
shear(:,1) = time;
shear(:,2) = CBFv.*pi.*(diam./4).^2;%Divide by 4 to convert to radius of 0.1
shear(:,2) = shear(:,2);%+0.6*rand(size(CBF(:,2)));
Press = 70.*ones(size(shear)); %Assume pressure is normal CPP = 70 mmHg

load("~/Documents/GitHub/CA_Hemodynamics/Code/Paramvals.mat")
%The paramvals used were: 
% paramvals = [0.1, 5.07, 2.5, 40, 25, 1, 1, 1, 0, 40, 0.5, 100000, 10.5]

%% 
% Define initial conditions
IC = [paramvals(2), paramvals(1), 0, 0, 0, 0, 0, 0, 0];

% Run model (note pressure is converted from mmHg to N/m^2
[t,y] = ode23(@(t, y) CereBRLSIM_FMD(t,y,paramvals,[time, Press/75], shear, [], []), time(4:end)', IC);


%scale diameters:
b1 = [ones(size(y(:,2))) y(:,2)]\diam(4:end);
if length(find(isnan(b1))>1)
    b1 = [ones(size(y(:,2))) y(:,2)]\diam(1:end-3);
end
sim_diam = [ones(size(y(:,2))) y(:,2)]*(b1);


fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot(t, sim_diam, 'k', 'linewidth',3), hold on, plot(time, diam, 'k:', 'linewidth',3) %don't have to convert because radius doesn't change with CVTR off
ylabel('Diameter (mm)')


legend('Scaled Diameter_{sim}','Diameter_{data}', 'FontSize',15);%, 'CBF_{data}')%, title(['Radius = ', num2str(paramvals(1))])

set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',25)
xlabel('Time (s)')

figure, plot(t, y(:,1))