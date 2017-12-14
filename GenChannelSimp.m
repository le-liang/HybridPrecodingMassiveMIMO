% Generate a simplified MU-MIMO channel based on uniform linear array (ULA)
% "The capacity optimality of beam steering in large mmWave MIMO channels"
% [H, Gain, At] = GenChannelSimp (Nt, K, Ncls)
% Input: Nt, number of TX antennas
%        K, number of single-antenna UEs
%        Ncls, number of clusters, single path per cluster
%        d, normalized antenna spacing
% Output: H, generated MU-MIMO channel
%         Gain, complex gain of each path
%         At, concatenation of array response vectors

% By Le Liang, UVic


function [H, Gain, At] = GenChannelSimp(Nt, K, Ncls, d)
if nargin <= 3
    d = 1/2;
end

% spread = pi/12;% [-15, 15]
spread = pi;
ang = (2*rand(Ncls, K)-1)*spread; % generated path angles, uniform distribution
Gain = (randn(Ncls, K) + j*randn(Ncls, K))/sqrt(2);

At = zeros(Nt, Ncls, K);
H = zeros(K, Nt);
for ik = 1 : K
    tmp = 0 : (Nt-1);
    tmp = tmp(:);
    kron1 = j*2*pi*d*tmp;
    kron2 = sin(ang(:, ik));
    kron2 = kron2.';
    kronr = kron(kron1, kron2);
    
    At(:, :, ik) = 1/sqrt(Nt) * exp(kronr);
    H(ik, :) = sqrt(Nt/Ncls) * Gain(:, ik).' * At(:, :, ik)';
end
    
    