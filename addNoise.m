function [vector] = addNoise(vector, mean, variance)
    noise = sqrt(variance) * randn(size(vector)) + mean;
    vector = vector + noise;
end