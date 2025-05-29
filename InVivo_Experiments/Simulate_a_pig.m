function [time, state, y] = Simulate_a_pig(CPP_Base, CPPchange)

%In this function we simulate the pig's reaction to change in CPP:

% Create CPP waveform over 2 hours:
t = [1:2*60*60]; %Time is in seconds

%First we make a logistic function from 1 to 1 which shapes the CPP change
logistic = 1./(1+exp(-1*sign(CPPchange)*0.002*(t-60*60)));

%Then we rescale to the acutual CPP values
C1 = CPP_Base+CPPchange/100*CPP_Base; %Gives the change
y = rescale(logistic,min(C1,CPP_Base),max(CPP_Base,C1)); %conversion from mmHg to N/cm^2 


y = y./75;
y = movmean(y, 2/0.01);

%----- Run the model! ----- %

%Load paramvals
paramvals = [15, 4, 2.5, 40, 25, 1, 0, 10, 0, 40, 0.5, 7000, 10.5, 0]
IC = [paramvals(2), paramvals(1), 0, 0, 0, 0, 0, 0, 0];


[time,state] = ode23(@(time, state) CereBRLSIM_Pigs(time,state,paramvals,[t', y']), t(4:end), IC);

end