%Simulate relationship between ETCO2 and CBF: 

%----- Run the model! ----- %
load("~/Documents/GitHub/CA_Hemodynamics/Code/Paramvals.mat")
%paramvals
paramvals = [0.10, 5.07, 1, 40, 1, 1, 10, 10, 0, 43, 0.5, 7000, 18, 0]

%CO2 - 
t = [0:0.1:15];
IC = [paramvals(2), paramvals(1), 0, 0, 0, 0, 0, 0, 0];


%P1: 
paramvals(10) = 46;
CO2up = linspace(paramvals(10),67,length(t));
P = linspace(89,101,length(t))+5*rand(length(t),1);
figure(1)
[time,state] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2up']), t(1:end), IC);
CO2down = linspace(paramvals(10),23,length(t));
P = linspace(89,65,length(t))+5*rand(length(t),1);
[time2,state2] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2down']), t(1:end), IC);
c = [0    0.4470    0.7410]
plot(CO2up, (state(:,1))/state(1,1) ,'*', 'color',c)
hold on, plot(CO2down, (state2(:,1))/state2(1,1), '*', 'color',c)

%P2
paramvals(10) = 42;
CO2up = linspace(paramvals(10),61,length(t));
P = linspace(77,92,length(t))+5*rand(length(t),1);
figure(1)
[time,state] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2up']), t(1:end), IC);
CO2down = linspace(paramvals(10),22,length(t));
P = linspace(77,57,length(t))+5*rand(length(t),1);
[time2,state2] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2down']), t(1:end), IC);
c = [0    0.4470    0.7410]
plot(CO2up, (state(:,1))/state(1,1) ,'*', 'color',c)
hold on, plot(CO2down, (state2(:,1))/state2(1,1), '*', 'color',c)



%P3
paramvals(10) = 41;
CO2up = linspace(paramvals(10),58,length(t));
P = linspace(90,142,length(t))+5*rand(length(t),1);
figure(1)
[time,state] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2up']), t(1:end), IC);
CO2down = linspace(paramvals(10),24,length(t));
P = linspace(90,63,length(t))+5*rand(length(t),1);
[time2,state2] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2down']), t(1:end), IC);
c = [0    0.4470    0.7410]
plot(CO2up, (state(:,1))/state(1,1) ,'*', 'color',c)
hold on, plot(CO2down, (state2(:,1))/state2(1,1), '*', 'color',c)

%P4
paramvals(10) = 42;
CO2up = linspace(paramvals(10),62,length(t));
P = linspace(89,101,length(t))+5*rand(length(t),1);
figure(1)
[time,state] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2up']), t(1:end), IC);
CO2down = linspace(paramvals(10),16,length(t));
P = linspace(89,73,length(t))+5*rand(length(t),1);
[time2,state2] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2down']), t(1:end), IC);
c = [0    0.4470    0.7410]
plot(CO2up, (state(:,1))/state(1,1) ,'*', 'color',c)
hold on, plot(CO2down, (state2(:,1))/state2(1,1), '*', 'color',c)


%P5
paramvals(10) = 43;
CO2up = linspace(paramvals(10),61,length(t));
P = linspace(86,103,length(t))+5*rand(length(t),1);
figure(1)
[time,state] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2up']), t(1:end), IC);
CO2down = linspace(paramvals(10),17,length(t));
P = linspace(86,65,length(t))+5*rand(length(t),1);
[time2,state2] = ode23(@(time, state) CereBRLSIM_CO2(time,state,paramvals,[t', P/75], [t', CO2down']), t(1:end), IC);
c = [0    0.4470    0.7410]
plot(CO2up, (state(:,1))/state(1,1) ,'*', 'color',c)
hold on, plot(CO2down, (state2(:,1))/state2(1,1), '*', 'color',c)


T = readmatrix('PETCO2_vs_CBFv.csv')
[sorted,sorti] = sort(T(:,1));
sorted2 = T(sorti,2);
hold on, plot(sorted, movmean(sorted2./100, 4), 'k', 'linewidth', 3)

xlabel('ETCO2 (mmHg)')
ylabel('Change in Flow')
set(gcf, 'color','white')
set(gca, 'box','off')
set(gca, 'fontsize',20)
