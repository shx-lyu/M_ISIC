function MMSE = MMSE_filters

MMSE.SMI_MMSE = @SMI_MMSE;
MMSE.Ideal_MMSE = @Ideal_MMSE;

end

% SMI_MMSE
function smi_rate = SMI_MMSE(s,Ry_hat_inv,Y,B,sourceNum,host_height,host_width)
% Sample-matrix-inversion MMSE
% Parameters: input: carriers:s; inverse of sample average autocorrelation matrix:Ry_hat_inv
%                    origin watermark:B; number of carriers:sourceNum
% Output: ber of extracted watermark:smi_rate
    b_hat = sign(s'*Ry_hat_inv*Y);
    ber_mat = (b_hat == B);
    smi_rate = length(find(ber_mat == 0))/(sourceNum*(host_height*host_width));
end

% Ideal_MMSE
function Ideal_rate = Ideal_MMSE(s,Ry_inv,Y,B,sourceNum,host_height,host_width)
% Minimum mean square error 
% Parameters: input: carriers:s; inverse of autocorrelation matrix:Ry_inv
%                    origin watermark:B; number of carriers:sourceNum
% Output: ber of extracted watermark:Ideal_rate
    b_hat = sign(s'*Ry_inv*Y);
    ber_mat = (b_hat == B);
    Ideal_rate = length(find(ber_mat == 0))/(sourceNum*host_height*host_width);
end

