function channelSamples = GenChannelSamples(nSamples, alpha, mu, ms, omegaS)
        omegaF = gamrnd(ms, omegaS, [1, nSamples]).^(alpha/2);
        channelSamples = gamrnd(mu, omegaF / mu).^(1/alpha);
end