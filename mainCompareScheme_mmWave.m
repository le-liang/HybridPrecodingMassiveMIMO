% DESCRIPTION
%   mmWave massive MU-MIMO under chain limitations, single-antenna UE
%   Hybrid precoding with ZF at baseband and Phase Reversal (PR) Preocidng at analog
%   Compare full-complexity ZF (FC-ZF), Hybrid, QHybrid, and B-MIMO numerically
%   Reproduce Fig. 3 in [1]

% REFERENCE
%   [1] L. Liang, W. Xu, and X. Dong, "Low-complexity hybrid precoding in massive 
%   mulituser MIMO systems," IEEE Wireless Communications Letters, vol. 3,
%   no. 6, pp. 653-656, Dec. 2014.

% By Le Liang, UVic, Mar. 25, 2014

tic; clear all; clc

Nt = 128;
K = 4; % UE number
Np = 10; % number of paths per user
B1 = 1; % quantized analog beamforming, up to B bits of precision
B2 = 2;

SNR = -30 : 5 : 0; 
nSNR = length(SNR);
channNum = 1e3;

rateZF = zeros(nSNR, 1); % FC-ZF
rateHyb = zeros(nSNR, 1);% W = ZF at baseband, F = PR at analog
rateHybQ1 = zeros(nSNR, 1);% Quantized hybrid precoding
rateHybQ2 = zeros(nSNR, 1);
rateBMIMO = zeros(nSNR, 1);% Multiuser beamspace MIMO precoder (BMIMO)

for isnr = 1 : nSNR
    
    P = 10^(SNR(isnr)/10);
    
    for ichannel = 1 : channNum
        
        [H, Gain, At] = GenChannelSimp(Nt, K, Np, 0.5); % mmWave channel
        % H = K x Nt, Gain = Np x K, At = Nt x Np x K
        
        %   =============  ZF preocidng, numerical  ================ 
        WtZF = H'*inv(H*H');
        WZF = WtZF*inv(sqrt(diag(diag(WtZF'*WtZF)))); % normalized columns
        rateZF(isnr) = rateZF(isnr) + CalRate(P/K*eye(K), H, WZF);
        
        %   ============ Hybrid precoding, numerical ===============
        for ik = 1 : K
            ph = - phase(H(ik,:));
            ph = ph(:);
            F(:,ik) =1/sqrt(Nt)* exp(j.*ph); % analog RF preprocessing
        end
        Fb = CalBDPrecoder(H*F);% digital baseband, same as inverse with column normalization
        wt = F*Fb;% aggregated precoder
        WPR = wt*inv(sqrt(diag(diag(wt'*wt))));   
        rateHyb(isnr) = rateHyb(isnr) + CalRate((P/K)*eye(K), H, WPR);% ZF-PRP
        
        %   ============= Quantized hybrid precoding ============
        FQPR1 = 1/sqrt(Nt) * Quant(B1, F);% analog RF
        wt = FQPR1*CalBDPrecoder(H*FQPR1);
        WQPR1 = wt*inv(sqrt(diag(diag(wt'*wt))));
        rateHybQ1(isnr) = rateHybQ1(isnr) + CalRate((P/K)*eye(K), H, WQPR1);
        
        WtQPR2 = 1/sqrt(Nt) * Quant(B2, F);
        wt = WtQPR2*CalBDPrecoder(H*WtQPR2);
        WQPR2 = wt*inv(sqrt(diag(diag(wt'*wt))));
        rateHybQ2(isnr) = rateHybQ2(isnr) + CalRate((P/K)*eye(K), H, WQPR2);

        % =========== Multiuser Beamspace MIMO precoder (B-MIMO) ============== %
        D = dftmtx(Nt);
        Hf = H*D;
        [maxVal, maxInd] = sort(diag(Hf'*Hf), 'descend'); % sort column with decreasing magnitude
        FBMIMO = D(:, maxInd(1:K));
        FbBMIMO = pinv(H*FBMIMO);
        Wt = FBMIMO*FbBMIMO;% analog x baseband precoding
        WBMIMO = Wt*inv(sqrt(diag(diag(Wt'*Wt))));
        rateBMIMO(isnr) = rateBMIMO(isnr) + CalRate((P/K)*eye(K), H, WBMIMO);
    end
    isnr   
end

rateZF = rateZF/channNum;
rateHyb = rateHyb/channNum;
rateHybQ1 = rateHybQ1/channNum;
rateHybQ2 = rateHybQ2/channNum;
rateBMIMO = rateBMIMO/channNum;

LineWidth = 1.5;
MarkerSize = 6;
figure
plot(SNR, abs(rateZF), 'k-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateHyb),'r-*', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateHybQ1), 'b-^', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateHybQ2), 'b-v', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(SNR, abs(rateBMIMO), 'm-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold off
legend('FC-ZF Precoding', 'Hybrid Precoding', 'Quantized Hybrid Precoding, B = 1',...
    'Quantized Hybrid Precoding, B = 2', 'B-MIMO Preocoding');
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bps/Hz)')
% title(sprintf('Nt = %d, K = %d,  Np = %d',Nt, K, Np))
grid

saveas(gcf, sprintf('MainCompareScheme-Nt%d-K%d-Np%d', Nt, K, Np)); % save current figure to file
toc


