% Quant(B, W) quantizes phases of each element in W up to B bits of precision
% By Le Liang, UVic, April 15, 2014

function r = Quant(B, W)
delta = 2*pi/2^B; % quantization interval
r = zeros(size(W, 1), size(W, 2));% ininitialize quantized matrix

for i1 = 1 : size(W, 1)
    for i2 = 1 : size(W, 2)
        ph = phase(W(i1, i2)); % ph in [-pi, pi]
        phq = floor(ph/delta)*delta +(mod(ph, delta) > delta/2)*delta ;% quantized phase
        r(i1, i2) = exp(j*phq);
    end
end
r = 1/sqrt(size(W, 1)) * r;
end