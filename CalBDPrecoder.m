%  BD precoder design subfunction r = CalBDPrecoder(G, K) --%
%  Input:  H, K*Nr x Nt,  aggregate channel, 
%          NsUE, # streams per user, K x 1
%  Output: r, Nt x K*Lr BD precoder

%  By Le Liang, UVic, Oct. 27, 2013

function r = CalBDPrecoder(H, NsUE) % no waterfilling
if nargin <= 1
    NsUE = ones(size(H,1)); % essentially ZF precoding with column normalization
end
[NR, Nt] = size(H);
K = length(NsUE);
Nr = NR/K;
F = [];% BD precoder initialization

%**** Standard BD design based on paper Q. H. Spencer, et al.,"Zero-forcing methods for downlink spatial multiplexing in multiuser MIMO channels," IEEE TSP, Feb. 2004. ***%
for iK = 1 : K
    Gj = H(((iK-1)*Nr+1):(iK*Nr),:);
    G_tilde = H;
    G_tilde(((iK-1)*Nr+1):(iK*Nr),:) = []; % Delete iK-user's channel, G_tilde
    [U, S, V] = svd(G_tilde);
    rv = rank(G_tilde);
    
    Htmp = Gj*V(:, (rv + 1):Nt );
    Lr = rank(Htmp);% # of data streams accommodated by iK user
    [Ut, St, Vt] = svd(Htmp);
    Ftmp = V( : , (rv + 1):Nt)*Vt(: , 1:NsUE(iK));
    F = [F Ftmp];
end
r = F;

