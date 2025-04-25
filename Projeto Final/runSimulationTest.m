close all; clear; clc;
numSamples = [5000];
numPoints = 51;
deltaAtt = [2.001];
alpha = [2];
mu = [1, 5];
snrdB = [-12];
numMontCarlo = 10^5;

runSimulationAndSaveDat('numSamples', numSamples, 'numPoints', numPoints, 'alpha', alpha, 'mu', mu, 'numMontCarlo', numMontCarlo, 'snrdB', snrdB, 'deltaAtt', deltaAtt);

