function JADE_ber = JADE(Y,m,B,host_height,host_width)
% JADE.Reference:Jean-Franccois Cardoso,et al."Blind beamforming for non Gaussian signals","High-order contrasts for independent component analysis"
    W = jadeR(Y,m);
    extractedWatermark = sign(W*Y);
    [JADE_ber,~] = BER(B,extractedWatermark,m,host_height,host_width);
end