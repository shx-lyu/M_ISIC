function smi_rate = SMI_MMSE(s,Ry_hat_inv,Y,B,sourceNum,host_height,host_width)
% Sample-matrix-inversion MMSE
% parameters: input: carriers:s; inverse of sample average autocorrelation matrix:Ry_hat_inv
%                    origin watermark:B; number of carriers:sourceNum
% output: ber of extracted watermark:smi_rate
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
    b_hat = sign(s'*Ry_hat_inv*Y);
    ber_mat = (b_hat == B);
    smi_rate = length(find(ber_mat == 0))/(sourceNum*(host_height*host_width));
end