function Ideal_rate = Ideal_MMSE(s,Ry_inv,Y,B,sourceNum,host_height,host_width)
% Minimum mean square error 
% parameters: input: carriers:s; inverse of autocorrelation matrix:Ry_inv
%                    origin watermark:B; number of carriers:sourceNum
% output: ber of extracted watermark:Ideal_rate
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
    b_hat = sign(s'*Ry_inv*Y);
    ber_mat = (b_hat == B);
    Ideal_rate = length(find(ber_mat == 0))/(sourceNum*host_height*host_width);
end