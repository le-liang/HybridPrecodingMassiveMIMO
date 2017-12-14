% DESCRIPTION
%   Massive MU-MIMO under chain limitations, single-antenna UE
%   Hybrid precoding with ZF at baseband and phase reversal at analog
%   Compare full-complexity ZF (FC-ZF) vs Hybrid
%   Reproduce Fig. 2 in [1]

% REFERENCE
%   [1] L. Liang, W. Xu, and X. Dong, "Low-complexity hybrid precoding in massive 
%   mulituser MIMO systems," IEEE Wireless Communications Letters, vol. 3,
%   no. 6, pp. 653-656, Dec. 2014.

% By Le Liang, UVic, Mar. 25, 2014

tic; clear; clc

Nt = 8;
K = 2; % UE number
B1 = 1; % quantized analog beamforming, up to B bits of precision
B2 = 2;

SNR = 0 : 5 : 30; 
nSNR = length(SNR);
channNum = 1e3;

rateZF = zeros(nSNR, 1); % FC-ZF
rateZFA = zeros(nSNR, 1);
rateHyb = zeros(nSNR, 1);% W = ZF at baseband, F = PR at analog
rateHybA = zeros(nSNR, 1);
rateQHyb1 = zeros(nSNR, 1);% Quantized PR, ZF at baseband
rateQhyb2 = zeros(nSNR, 1);

for isnr = 1 : nSNR
    
    P = 10^(SNR(isnr)/10);
    
    % ====================================================================
    % ===================== Analytical result ============================
    % ====================================================================
    % ==========  ZF analytical results =============
    s = 0;
    for k = 1 : (Nt)
        rho = P/K;
        s = s + exp(1/rho)*log2(exp(1)) * expintn_noMaple(1/rho, k); % use maple for expintn()
        % s = s + exp(1/rho)*log2(exp(1)) * expintn_noMaple(1/rho, k); % if no maple installed
    end    
    rateZFA(isnr) = K*s;
    
    % ========= Hybrid precoding analytical results  ==========
    rateHybA(isnr) = K*log2(1 + (pi/4) * (P*Nt)/K );
    
    % ====================================================================
    % ===================== Simulation result ============================
    % ====================================================================   
    for ichannel = 1 : channNum
        
        H = (randn(K, Nt) + j*randn(K, Nt))/sqrt(2);
        
        %   =============  ZF preocidng, numerical  ================ 
        WtZF = H'*inv(H*H');
        WZF = WtZF*inv(sqrt(diag(diag(WtZF'*WtZF)))); % normalized columns
        rateZF(isnr) = rateZF(isnr) + CalRate(P/K*eye(K), H, WZF);
        
        %   ============ Hybrid, numerical ===============
        F = 1/sqrt(Nt)*exp(j.*angle(H))';
        Fb = CalBDPrecoder(H*F);% baseband, same as inverse with column normalization
        wt = F*Fb;% aggregate precoder
        WPR = wt*inv(sqrt(diag(diag(wt'*wt))));   
        rateHyb(isnr) = rateHyb(isnr) + CalRate((P/K)*eye(K), H, WPR);% ZF-PRP
        
        %   ============= Quantized ZF-PR precoding ============
        FQPR1 = Quant(B1, F);% analog RF
        wt = FQPR1*CalBDPrecoder(H*FQPR1);
        WQPR1 = wt*inv(sqrt(diag(diag(wt'*wt))));
        rateQHyb1(isnr) = rateQHyb1(isnr) + CalRate((P/K)*eye(K), H, WQPR1);
        
        WtQPR2 = Quant(B2, F);
        wt = WtQPR2*CalBDPrecoder(H*WtQPR2);
        WQPR2 = wt*inv(sqrt(diag(diag(wt'*wt))));
        rateQhyb2(isnr) = rateQhyb2(isnr) + CalRate((P/K)*eye(K), H, WQPR2);      
    end
    isnr   
end

rateZF = rateZF/channNum;
rateHyb = rateHyb/channNum;
rateQHyb1 = rateQHyb1/channNum;
rateQhyb2 = rateQhyb2/channNum;

LineWidth = 1.5;
MarkerSize = 6;
figure
plot(SNR, abs(rateZF), 'k-x', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateZFA), 'bo', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);
hold on
plot(SNR, abs(rateHyb),'r-*', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateHybA), 'gd', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateQHyb1), 'b-^', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateQhyb2), 'b-v', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold off
legend('FC-ZF, simulation', 'FC-ZF, analytical', ...
    'PZF, simulation', 'PZF, analytical', ...
    'Quantized PZF, B = 1','Quantized PZF, B = 2');
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bps/Hz)')
title(sprintf('Nt = %d, K = %d', Nt, K))
grid
% saveas(gcf, sprintf('MassiveCompareScheme-Nt%d-K%d', Nt, K)); 
% saveas(gcf, sprintf('MassiveCompareScheme-Nt%d-K%d-PhaseOfZf', Nt, K));

toc


