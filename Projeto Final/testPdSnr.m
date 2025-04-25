%% params
close all; clear; clc;
% warning('off','all')
numSamples = [1000];
deltaAtt = [2];
alpha = [2];
mu = [2];
ms = [2 ];
snrdB = [-21:3:9];
numMontCarlo = 10^5;
isMobile = [true];
numBins=40;
rm=100;
d0=1;
pFaStar = 0.1;
pDStar = 0.91;
maxIter = 10;

sampleFactor = 100;

maxDist = (d0 + rm) / d0;


delta = 2 * pi / numBins;


simPoints = length(snrdB);
pdSim = zeros(1, simPoints);
pfSim = zeros(1, simPoints);
pdTheo = zeros(1, simPoints);
pdAsymp = zeros(1, simPoints);
nSamples = zeros(1, simPoints);
errors = zeros(1, simPoints);

%[averageTheo, averageSim] = simulator.ReturnAveragePower();

%% Find N Samples
for ii = 1:simPoints
    pTx = 10^(snrdB(ii) / 10) * maxDist^(deltaAtt);

    simulator = PhaseDiffSimulatorV3(numMontCarlo, ...
        50, ...
        pTx, ...
        0, ...
        numBins, ...
        sampleFactor, ...
        alpha+randn(1)*1e-2, ...
        mu+randn(1)*1e-2, ...
        ms+randn(1)*1e-2, ...
        deltaAtt+randn(1)*1e-2, ...
        rm, ...
        d0, ...
        isMobile);

    % Initialization
    iter = 0;
    maxIter = 100;  % Set a reasonable max iteration limit
    err = inf;
    tolerance = 0.01;  % Convergence criterion
    pD = pDStar;
    disp('Iteration | numSamples | pdHat   | pD     | err'); % Header for logging
    disp('-------------------------------------------------');
    while abs(err)/pDStar >= tolerance && iter < maxIter
        % Compute the number of samples without correction
        numSamples = ceil((simulator.estimateNumberSamples(pFaStar, pD)));
        numSamples = min(1e6, numSamples);
        % Compute the threshold
        threshold = sqrt((2 * pi - delta) / (4 * pi * numSamples)) * qfuncinv(pFaStar);

        % Update simulator parameters
        simulator.threshold = threshold;
        simulator.numSamples = numSamples;

        % Compute the estimated Pd
        pdHat = real(simulator.computeTheoreticalPD());%

        % Compute the error
        err = pDStar - pdHat;

        % Display iteration log
        fprintf('%4d      | %8d  | %.4f  | %.4f  | %.4f\n', iter, numSamples, pdHat, pD, err);
        % Update pD cautiously to avoid overshooting
        if abs(err) > tolerance * pDStar
            pD = pD + 500 * err;  % Apply partial correction
        end

        iter = iter + 1;  % Proper iteration update
    end

    % Final simulation run
    practicalSamples = max(2*numBins, numSamples);
    threshold = sqrt((2 * pi - delta) / (4 * pi * practicalSamples)) * qfuncinv(pFaStar);

    % Update simulator parameters
    simulator.threshold = threshold;
    simulator.numSamples = practicalSamples;
    nSamples(ii) = practicalSamples;
    [pdSim(ii), pfSim(ii), ~, ~] = simulator.runSimulation();
    pdTheo(ii) = real(simulator.computeTheoreticalPD());
    pdAsymp(ii) = real(simulator.computeAsympPD(pFaStar));

end



% Plotting
figure(1)
hold on;
plot(snrdB, pdSim, 'DisplayName', ...
    'simulated-values', ...
    'Marker', 'square', 'LineWidth', 2, 'LineStyle', '--');
xlabel('SNR');
ylabel('P_{D}');
hold on;
grid on;

plot(snrdB, pdTheo, 'DisplayName',...
    'theoretical_values', ...
    'LineWidth', 2);

plot(snrdB, pdAsymp, 'DisplayName',...
    'Asymptotic_values', ...
    'LineWidth', 2);

axis([snrdB(1)  snrdB(length(snrdB)) 0 1]);
legend('show', 'Location','southwest')
hold off;

%% Fixed N Samples
for ii = 1:simPoints
    pTx = 10^(snrdB(ii) / 10) * maxDist^(deltaAtt);
    % pTx = 10^(snrdB(ii) / 10);

    simulator = PhaseDiffSimulatorV3(numMontCarlo, ...
        50, ...
        pTx, ...
        0, ...
        numBins, ...
        sampleFactor, ...
        alpha+randn(1)*1e-2, ...
        mu+randn(1)*1e-2, ...
        ms+randn(1)*1e-2, ...
        deltaAtt+randn(1)*1e-2, ...
        rm, ...
        d0, ...
        isMobile);

    % Final simulation run
    practicalSamples = 2000;
    threshold = sqrt((2 * pi - delta) / (4 * pi * practicalSamples)) * qfuncinv(pFaStar);

    % Update simulator parameters
    simulator.threshold = threshold;
    simulator.numSamples = practicalSamples;

    [pdSim(ii), pfSim(ii), ~, ~] = simulator.runSimulation();
    pdTheo(ii) = real(simulator.computeTheoreticalPD());
    pdAsymp(ii) = real(simulator.computeAsympPD(pFaStar));

end


% Plotting
figure(1)
hold on;
plot(snrdB, pdSim, 'DisplayName', ...
    'simulated-values', ...
    'Marker', 'square', 'LineWidth', 2, 'LineStyle', '--');
xlabel('SNR');
ylabel('P_{D}');
hold on;
grid on;

plot(snrdB, pdTheo, 'DisplayName',...
    'theoretical_values', ...
    'LineWidth', 2);

plot(snrdB, pdAsymp, 'DisplayName',...
    'Asymptotic_values', ...
    'LineWidth', 2);

axis([snrdB(1)  snrdB(length(snrdB)) 0 1]);
legend('show', 'Location','southwest')
hold off;

