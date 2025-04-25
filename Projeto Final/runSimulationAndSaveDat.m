function runSimulationAndSaveDat(varargin)

% Parse input arguments
p = inputParser;
addParameter(p, 'numMontCarlo', 10^5);
addParameter(p, 'numSamples', 5000);
addParameter(p, 'numBins', 40);
addParameter(p, 'numPoints', 31);
addParameter(p, 'snrdB', -12);
addParameter(p, 'alpha', 0.501);
addParameter(p, 'mu', 2.01);
addParameter(p, 'deltaAtt', 2);
addParameter(p, 'd0', 1);
addParameter(p, 'rm', 50);

parse(p, varargin{:});

numMontCarlo = p.Results.numMontCarlo;
numSamples = p.Results.numSamples;
numBins = p.Results.numBins;
numPoints = p.Results.numPoints;
snrdB = p.Results.snrdB;
alpha = p.Results.alpha;
mu = p.Results.mu;
deltaAtt = p.Results.deltaAtt;
d0 = p.Results.d0;
rm = p.Results.rm;

if isnumeric(numMontCarlo) && isscalar(numMontCarlo)
    numMontCarlo = numMontCarlo(:)';
end
if isnumeric(numSamples) && isscalar(numSamples)
    numSamples = numSamples(:)';
end
if isnumeric(numBins) && isscalar(numBins)
    numBins = numBins(:)';
end
if isnumeric(numPoints) && isscalar(numPoints)
    numPoints = numPoints(:)';
end
if isnumeric(snrdB) && isscalar(snrdB)
    snrdB = snrdB(:)';
end
if isnumeric(alpha) && isscalar(alpha)
    alpha = alpha(:)';
end
if isnumeric(mu) && isscalar(mu)
    mu = mu(:)';
end
if isnumeric(deltaAtt) && isscalar(deltaAtt)
    deltaAtt = deltaAtt(:)';
end
if isnumeric(d0) && isscalar(d0)
    d0 = d0(:)';
end
if isnumeric(rm) && isscalar(rm)
    rm = rm(:)';
end

f = waitbar(0, '0 %', 'Name', 'running simulation...');
totalSim = length(numMontCarlo) * length(numSamples) * ...
    length(numBins) * length(numPoints) * length(snrdB) * ...
    length(alpha) * length(mu) * length(deltaAtt) * ...
    length(d0) * length(rm);

currentSim = 0;
colorIndex = 1;

currentDateFolder = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
mkdir(currentDateFolder);

colors = colormap(parula(totalSim));
% not indented for readability reasons
for i = 1:length(numMontCarlo)
for j = 1:length(numSamples)
for k = 1:length(numBins)
for l = 1:length(numPoints)
for m = 1:length(snrdB)
for n = 1:length(alpha)
for o = 1:length(mu)
for q = 1:length(deltaAtt)
for r = 1:length(d0)
for s = 1:length(rm)
    
    runSingleSimulationAndSaveDat(numMontCarlo(i), numSamples(j), numBins(k), numPoints(l), ...
    snrdB(m), alpha(n), mu(o), deltaAtt(q), d0(r), rm(s), colors(colorIndex, :));
    
    colorIndex = colorIndex + 1;

end; end; end; end; end; end; end; end; end; end
delete(f);

    function runSingleSimulationAndSaveDat(numMontCarlo, ... 
                                           numSamples, ...   
                                           numBins, ...  
                                           numPoints, ...    
                                           snrdB, ...    
                                           alpha, ...    
                                           mu, ...   
                                           deltaAtt, ... 
                                           d0, ...   
                                           rm, ...   
                                           color)
        
        sampleFactor = 100;
        
        maxDist = (d0 + rm) / d0;
        pTx = 10^(snrdB / 10) * maxDist^(deltaAtt);
        
        delta = 2 * pi / numBins;
        pfaVec = logspace(-6, -1e-6, numPoints);
        thresholdVec = sqrt((2 * pi - delta) / (4 * pi * numSamples)) ...   
            .* qfuncinv(pfaVec);
        simulator = PhaseDiffSimulator(numMontCarlo, ...
                                         numSamples, ...
                                         pTx, ...
                                         0, ...
                                         numBins, ...
                                         sampleFactor, ...
                                         alpha, ...
                                         mu, ...
                                         deltaAtt, ...
                                         rm, ...
                                         d0);
        simPoints = length(thresholdVec);
        pdSim = zeros(1, simPoints);
        pfSim = zeros(1, simPoints);
        pdTheo = zeros(1, simPoints);
        
        for ii = 1:simPoints
            currentSim = currentSim + 1;
            threshold = thresholdVec(ii);
            
            simulator.threshold = threshold;
            [pdSim(ii), pfSim(ii), ~, ~] = simulator.runSimulation();
            pdTheo(ii) = real(simulator.computeTheoreticalPD());
            prct = currentSim / (totalSim * simPoints);
            waitbar(prct,f,sprintf('%.2f %%',prct*100))
        end
        
        
        % Create file names within the folder
        simulatedFilename = fullfile(currentDateFolder,...
            sprintf('simulated_values_%d_%d_%d_%0.1f_%0.1f_%d_%d_%d.dat', ...
            numSamples, numBins, snrdB, alpha, mu, deltaAtt, d0, rm));
        theoreticalFilename = fullfile(currentDateFolder,...    
            sprintf('theoretical_values_%d_%d_%d_%0.1f_%0.1f_%d_%d_%d.dat', ...
            numSamples, numBins, snrdB, alpha, mu, deltaAtt, d0, rm));
        aucFilename = fullfile(currentDateFolder,...    
            sprintf('auc_values_%d_%d_%d_%0.1f_%0.1f_%d_%d_%d.dat', ...
            numSamples, numBins, snrdB, alpha, mu, deltaAtt, d0, rm));
        
        aucSim = trapz(pfSim, pdSim);
        aucTheo = trapz(pfaVec, pdTheo);

        % Save simulated and theoretical values in .dat files
        saveDat(simulatedFilename, pfSim, pdSim);
        saveDat(theoreticalFilename, pfaVec, pdTheo);
        saveDat(aucFilename, aucSim, aucTheo);
        
        % Plotting
        figure(1)
        hold on;
        plot([0, pfSim], [0, pdSim], 'DisplayName', ... 
            sprintf('simulated_values_%d_%d_%d_%0.1f_%0.1f_%d_%d_%d', ...
            numSamples, numBins, snrdB, alpha, mu, deltaAtt, d0, rm), ...
            'Marker', 'square', 'LineWidth', 2, 'LineStyle', '--', 'Color', color);
        xlabel('P_{FA}');
        ylabel('P_{D}');
        hold on;
        grid on;
        
        plot([0, pfaVec], [0, pdTheo], 'DisplayName',...    
            sprintf('theoretical_values_%d_%d_%d_%0.1f_%0.1f_%d_%d_%d', ...
            numSamples, numBins, snrdB, alpha, mu, deltaAtt, d0, rm),...
            'LineWidth', 2, 'Color', color);

        axis([0 1 0 1]);
        legend('show', 'Location','southwest')
        hold off;
    end

    function saveDat(filename, x, y)
        fid = fopen(filename, 'w');
        fprintf(fid, '%f %f\n', [x; y]);
        fclose(fid);
    end
end