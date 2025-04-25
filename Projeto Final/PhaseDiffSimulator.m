classdef PhaseDiffSimulator
    properties
        numMonteCarlo
        numSamples
        pTx
        threshold
        numBins
        sampleFactor
        thetaDemod
        alpha
        mu
        attenuationCoeff
        rm
        d0
        channelGeneratorObj
    end

    methods
        function obj = PhaseDiffSimulator(numMonteCarlo,...
                numSamples,...
                pTx,...
                threshold,...
                numBins,...
                sampleFactor,...
                alpha,...
                mu,...
                attenuationCoeff,...
                rm,...
                d0)
            obj.numMonteCarlo = numMonteCarlo;
            obj.numSamples = numSamples;
            obj.pTx = pTx;
            obj.threshold = threshold;
            obj.numBins = numBins;
            obj.sampleFactor = sampleFactor;
            obj.alpha = alpha;
            obj.mu = mu;
            obj.attenuationCoeff = attenuationCoeff;
            obj.rm = rm;
            obj.d0 = d0;
            obj.channelGeneratorObj = ChannelModel(alpha, mu, ...
                attenuationCoeff, rm, d0, pTx);
        end

        function [estPD, estPF, testH0, testH1] = runSimulation(obj)

            numSymbols = round(obj.numSamples / obj.sampleFactor);

            sampleRatio = 1;

            sampleFrequency = obj.sampleFactor * sampleRatio;  % Use obj.sampleFactor

            oscillatorFrequency = sampleFrequency / 4;

            % useful constants
            contPD = 0;
            contPF = 0;
            testH0 = zeros(1, obj.numMonteCarlo/2);  % Use obj.numMonteCarlo
            testH1 = zeros(1, obj.numMonteCarlo/2);  % Use obj.numMonteCarlo

            dataLength = obj.sampleFactor * numSymbols;  % Use obj.sampleFactor and obj.numSymbols
            tMod = 0: 1/sampleFrequency : (dataLength - 1)*1/sampleFrequency;
            modBaseFunctions = exp(1i * 2 * pi * oscillatorFrequency .* tMod);
            obj.thetaDemod = 2 * pi * oscillatorFrequency / sampleFrequency;  % Use sampleFrequency
            % sqrtSignalPowerLin = sqrt(obj.pTx);
            % sqrtSignalPowerLin = sqrt(1);
            % First Loop Tx is not enabled
            parfor i = 1 : obj.numMonteCarlo / 2
                noise = 1/sqrt(2) * randn(1, dataLength) ...
                    + 1i * 1/sqrt(2) * randn(1, dataLength);
                receivedData  =  noise;
                testVar = obj.performStatisticalTest(receivedData);
                testH0(i) = testVar;
                if testVar >= obj.threshold
                    contPF = contPF + 1;
                end
            end

            parfor i = 1 : obj.numMonteCarlo / 2
                data = obj.generateBPSKSymbols(numSymbols);
                sampledData = zeros(1, dataLength);
                for j = 1 : obj.sampleFactor
                    sampledData(j : obj.sampleFactor : end) = data;
                end
                modulatedData = sampledData .* modBaseFunctions;
                signalPower = sum(real(modulatedData).^2 + imag(modulatedData).^2) / dataLength;
                modulatedData = modulatedData / sqrt(signalPower) * sqrt(1);
                noise = 1/sqrt(2) * randn(1, dataLength) ...
                    + 1i * 1/sqrt(2) * randn(1, dataLength);
                x = obj.channelGeneratorObj.ReturnChannelSamples(1);
                receivedData  = modulatedData * x + noise;
                testVar = obj.performStatisticalTest(receivedData);
                testH1(i) = testVar;
                if testVar >= obj.threshold
                    contPD = contPD + 1;
                end
            end

            estPF = contPF / (obj.numMonteCarlo/2);
            estPD = contPD / (obj.numMonteCarlo/2);
        end



        function pd = computeTheoreticalPD(obj)

            [cn, Cn, cp, Cp, dm, Dm, dq, Dq, kappa, z] =...
                obj.channelGeneratorObj.ReturnChannelFoxHParams();


            delta = 2 * pi / obj.numBins;

            a = pi*sqrt(obj.numSamples)/(2*sqrt(2-delta/pi));

            b = 2*obj.threshold*sqrt(obj.numSamples)/sqrt(2-delta/pi);

            x = z.*sqrt(2)./a;

            y = -b^2/2;

            numTerms = length(z);
            auxSum = zeros(1, numTerms);
            for i = 1: numTerms
                if isempty(cp)
                    cp2 = [];
                    Cp2 = [];
                else
                    cp2 = cp(i, :);
                    Cp2 = Cp(i, :);
                end
                cn2 = cn(i, :);
                Cn2 = Cn(i, :);
                dm2 = dm(i, :);
                Dm2 = Dm(i, :);
                dq2 = dq(i, :);
                Dq2 = Dq(i, :);


                h1 = obj.computeFirstBiFox(x(i), y, cn2, Cn2, cp2, Cp2, dm2, ...
                    Dm2, dq2, Dq2);


                h2 = obj.computeSecodBiFox(x(i), y, cn2, Cn2, cp2, Cp2, dm2, ...
                    Dm2, dq2, Dq2);

                auxSum(i) = kappa(i)*(b*sqrt(2)/2*h1 + h2);

            end
            pd = 1 - exp(-b^2/2)/2*sum(auxSum);

        end

    end

    methods (Access = private)
        function testVar = performStatisticalTest(obj, receivedData)
            signalPhase = atan2(imag(receivedData), real(receivedData));

            phaseDiff = mod((signalPhase(2:end) - signalPhase(1:end-1)), 2 * pi);

            [xEmpirical, yEmpirical] = obj.estimateNumericalPDF(phaseDiff, obj.numBins);  % Use obj.numBins

            testVar = sum(yEmpirical .* cos(obj.thetaDemod - xEmpirical) * 2 * pi / obj.numBins);  % Use obj.thetaDemod and obj.numBins
        end

        function symbols = generateBPSKSymbols(~, num)
            normConstellation = [1 -1];
            symbols = normConstellation(randi([1, length(normConstellation)], 1, num));
        end

        function [xo, yo] = estimateNumericalPDF(~, samples, bins)
            [n, xout] = hist(samples, bins);
            n = n ./ sum(n) ./ (xout(2) - xout(1));
            xo = xout;
            yo = n;
        end

        function [result] = computeFirstBiFox(~, x, y, cn2, Cn2, cp2, Cp2, dm2, ...
                Dm2, dq2, Dq2)

            an1 = 0; alphan1 = 1/2; An1 = 1;
            ap1 = []; alphap1 = []; Ap1 = [];
            bq1 = []; betaq1 = []; Bq1 = [];

            en3 = []; En3 = []; ep3 = []; Ep3 = [];
            fm3 = 0; Fm3 = 1; fq3 = -1/2; Fq3 = 1;

            result = BiFoxH(an1, alphan1, An1, ap1, alphap1, Ap1,...
                bq1, betaq1 ,Bq1, ...
                cn2, Cn2, cp2, Cp2, ...
                dm2, Dm2, dq2, Dq2, ...
                en3, En3, ep3, Ep3, ...
                fm3, Fm3, fq3, Fq3, ...
                x, y);
        end

        function [result] = computeSecodBiFox(~, x, y, cn2, Cn2, cp2, Cp2, dm2, ...
                Dm2, dq2, Dq2)

            an1 = 1/2; alphan1 = 1/2; An1 = 1;
            ap1 = []; alphap1 = []; Ap1 = [];
            bq1 = []; betaq1 = []; Bq1 = [];

            en3 = []; En3 = []; ep3 = []; Ep3 = [];
            fm3 = 0; Fm3 = 1; fq3 = 1/2; Fq3 = 1;

            result = BiFoxH(an1, alphan1, An1, ap1, alphap1, Ap1,...
                bq1, betaq1 ,Bq1, ...
                cn2, Cn2, cp2, Cp2, ...
                dm2, Dm2, dq2, Dq2, ...
                en3, En3, ep3, Ep3, ...
                fm3, Fm3, fq3, Fq3, ...
                x, y);
        end

    end
end
