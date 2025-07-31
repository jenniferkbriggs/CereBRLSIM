%Run Experiment to model dCA
clear all
close all
clc

lagon = 1; %whether or not to lag between ABP to match blood velocity
filename = '07'
Data = readmatrix(['SQ_ABP.csv']); %First row is time, second row is ABP in mmHg
Flow = readmatrix(['SQ_v.csv']);   %First row is time, second row is blood velocity

%-- set paramvals:
paramvals = [0.10, 5.07, 2.5, 40, 25, 1, 1, 1, 0, 40, 0.5, 7000, 12.5, 0];

%interpolate data
Datasmoothed = movmean(Data(:,2),5);              
t = [Data(1,1):0.01:Data(end,1)];
[datasorted, indx] = sort(Data(:,1), 'descend');
if length(unique(datasorted))<length(Data(:,1))
    indx2 = find(diff(datasorted)==0)
    datasorted(indx2) = [];
    indx(indx2) = [];
end
y = interp1(datasorted, Datasmoothed(indx), t);
y = movmean(y, 6);
paramvals(14) = 1;


%interpolate flow:
Datasmoothed = movmean(Flow(:,2), 5);


%make sure all data are ordered
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
Flow = interp1(Flow(:,1),Flow(:,2), t);

%remove nans
bad = find(isnan(Flow));
t(bad) = [];
Flow(bad) = [];
y(bad) = [];



% Pressure is measured at the finger so the flow and pressure may be lagged
% from eachother. Find the max cross correlation and then shift by that
% lag. 
res = mean(diff(t))
[x,lag] = xcov(Flow, y, round(80/res)); %max lag is 10 seconds
[mx, indx] = max(x);
lag_s = lag(indx)*res; %this is the lag in seconds: 
% 
% % %lag flow backwards:
if lagon
t(1:-lag(indx)) = [];
y(1:-lag(indx)) = [];

Flow(end+lag(indx)+1:end) = [];
end


%remove nans
nanindx = find(isnan(Flow));
t(nanindx) = [];
y(nanindx) = [];
Flow(nanindx) = [];


IC = [paramvals(2), paramvals(1), 0, 0, 0, 0, 0, 0, 0];

[time,state] = ode23(@(time, state) CereBRLSIM_dCA(time,state,paramvals,[t', y'/75]), t(4:end), IC);

x = 0.1
rad_perc = (squeeze(state(:,2)))./squeeze(state(1,2)); %percent change of the radius
rad_perc = rad_perc-1; %quantify if it is getting bigger or smaller.  - this is what we will rescale for the MCA 
MCA_rad = (rad_perc.*x+1);
sim_cbfv = squeeze(state(:,1,:))./(MCA_rad.^2*pi);


%Scale
mdl = fitlm((sim_cbfv), (Flow(4:end)));



%plot model vs data
fig = figure, 
fig.Position = [-209 1482 370 364];
fig.Units = 'pixels'
plot((time-time(1)), mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*sim_cbfv, 'k:', 'linewidth',3), hold on,
plot((t-t(1)), Flow, 'k', 'linewidth',3)
legend('Scaled CBv_{sim}', 'CBv_{data}')%, title(['Radius = ', num2str(paramvals(1))])
set(gca, 'box','off')
set(gcf,'color','white')
set(gca, 'fontsize',15)
xlabel('Time (s)')
ylabel('CBv (cm/s)')

