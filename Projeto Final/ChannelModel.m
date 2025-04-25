classdef ChannelModel
    properties
        alpha
        mu
        delta
        rm
        d0
        maxDist
        pTx
    end

    methods
        function obj = ChannelModel(alpha, mu, delta, rm, d0, pTx)

            obj.alpha = alpha;
            obj.mu = mu;
            obj.delta = delta;
            obj.rm = rm;
            obj.d0 = d0;
            obj.maxDist = (rm + d0)/d0;
            obj.pTx = pTx;

        end

        function [an, An, ap, Ap, bm, Bm, bq, Bq, kappa, z] = ReturnChannelFoxHParams(obj)

            an = [1];
            An = [1];
            ap = [];
            Ap = [];
            bm = [obj.mu];
            Bm = [2/obj.alpha];
            bq = [0];
            Bq = [1];

            omega = obj.pTx .* obj.maxDist^(-obj.delta);
            z = obj.mu^(2/obj.alpha)/omega;

            kappa = 1 / (gamma(obj.mu));

        end

        function samples = ReturnChannelSamples(obj, numSamples)
            omega = obj.pTx*obj.maxDist.^(-obj.delta);
            samples = obj.generateChannelSamples(omega, numSamples);

        end

    end

    methods (Access = private)

        function channelSamples = generateChannelSamples(obj, omega, numSamples)
            omegaF = omega.^(obj.alpha/2);
            channelSamples = gamrnd(obj.mu, omegaF / obj.mu, [1, numSamples]).^(1/obj.alpha);
        end

    end
end