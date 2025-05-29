% Reproducing KleinDepreiteri2022 Quadraphic Curve study
clear all
clc

% Set random seed for reproducibility
rng(1)

% Number of pigs
numpigs = 200;


% Initialize baseline CPP
CPP_Base = normrnd(68.2, 10, [1, numpigs]); % Table 1 (approximate range of CPP increasing range from paper)

% Initialize variables
Flow = nan(numpigs, 2);
Radius = nan(numpigs, 2);
Pressure = nan(numpigs, 2);


% Select hypotensive pigs
hypo_pig = [1:numpigs/2];

% Select hypertensive pigs
hyper_pig = [numpigs/2+1:numpigs];

% Precompute CPP changes for hypotensive and hypertensive pigs

%decrease average percent change to fill in the center of the curve better
CPPchange_hypo = normrnd(-84.6, 25, [1, length(hypo_pig)]); % Randomly calculate the percent change of CPP
CPPchange_hypo(abs(CPPchange_hypo) > 100) = 99.;

CPPchange_hyper = normrnd(95.3, 25, [1, length(hyper_pig)]); % Randomly calculate the percent change of CPP


% Use a temporary struct array to store results within the parfor loop
results(numpigs) = struct('Flow', [], 'Radius', [], 'Pressure', [], 'Flow_all', [], 'Rad_all', [], 'Myo_all', [], 'Press_all', []);

for i = hypo_pig
    % Determine if the pig is hypotensive or hypertensive
        % Hypotensive pig
    CPPchange_idx = find(hypo_pig == i, 1);
    CPPchange = CPPchange_hypo(CPPchange_idx);
   

    [time, state, y] = Simulate_a_pig(CPP_Base(i), CPPchange);
    results(i).Flow = state([1, end], 1);
    results(i).Radius = state([1, end], 2);
    results(i).Pressure = y([1, end]) * 75.0062;

    results(i).Flow_all = state(:, 1);
    results(i).Rad_all = state(:, 2);
    results(i).Myo_all = state(:, 3);
    results(i).Press_all = y;
end
disp('Running hypertensive pigs')

for i = hyper_pig
        % Hypertensive pig
    CPPchange_idx = find(hyper_pig == i, 1);
    CPPchange = CPPchange_hyper(CPPchange_idx);

    [time, state, y] = Simulate_a_pig(CPP_Base(i), CPPchange);
    results(i).Flow = state([1, end], 1);
    results(i).Radius = state([1, end], 2);
    results(i).Pressure = y([1, end]) * 75.0062;

    results(i).Flow_all = state(:, 1);
    results(i).Rad_all = state(:, 2);
    results(i).Myo_all = state(:, 3);
    results(i).Press_all = y;
end

disp('Running hypotensive pigs')

% Transfer results from struct array to original arrays
for i = 1:numpigs
    Flow(i, :) = results(i).Flow;
    Radius(i, :) = results(i).Radius;
    Pressure(i, :) = results(i).Pressure;

    Flow_all(i, 1:length(results(i).Flow_all)) = results(i).Flow_all';
    Rad_all(i, 1:length(results(i).Rad_all)) = results(i).Rad_all';
    Myo_all(i, 1:length(results(i).Myo_all)) = results(i).Myo_all';
    Press_all(i, 1:length(results(i).Press_all)) = results(i).Press_all';
end



try % if data are present, plot and load them. 


% Load data
RBCFlux = readmatrix('RBCFlux_medium.csv');
Diameter = readmatrix('Diameter_medium.csv');

% Plot results
figure, plot(Pressure(:, 2), (Flow(:, 2) - Flow(:, 1)) ./ Flow(:, 1), '*')
ylabel([{'Change in Flow'}; {'% \Delta baseline'}]), xlabel('CPP')
hold on, plot(RBCFlux(:, 1), RBCFlux(:, 2), 'k', 'linewidth', 3)
legend('Flow Simulation', 'RBC Flux Data')
set(gca, 'FontSize', 15)
set(gcf, 'Color', 'white')
set(gca, 'box', 'off')

end