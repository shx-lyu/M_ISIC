function [BERs,iterate_times] = M_IGLS(Y,Ry_hat_inv,B,P,threshold,sourceNum,host_height,host_width,Bc_hat)
% Do P times M_IGLS,compute P bers
% parameters: input: embedded matrix:Y; prewhitening matrix:Ry_hat_inv;
%                    origin watermark:B; reinitialization times:P;
%                    threshold of iterating times:threshold;
%                    number of carriers:sourceNum
%                    pre-initialization:Bc_hat
% output: 1xP vector of BER:BERs; M_IGLS's actual iterating times:iterate_times
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
    k = 1;
    BERs = zeros(1,P);
    iterate_times = zeros(1,P);
    while(k<=P)
        [extractedWatermark,iterate_times(k)] = extract(Y,Ry_hat_inv,threshold,Bc_hat{1,k});
        [BERs(k),~] = BER(B,extractedWatermark,sourceNum,host_height,host_width);
        k=k+1;
    end
end

function [extractedWatermark,d] = extract(Y,Ry_hat_inv,threshold,B_hat)
% M_IGLS.Reference:M. Li, M. K. Kulhandjian, D. A. Pados, S. N. Batalama, and M. J.Medley, ¡°Extracting spread-spectrum hidden data from digital media,¡±IEEE transactions on information forensics and security, vol. 8, no. 7,pp. 1201¨C1210, 2013.
    d = 0;
    while(1)
        d = d+1;
        V_hat = (Y*B_hat')/(B_hat*B_hat');
        B_now = sign(((V_hat'*Ry_hat_inv*V_hat)+(1e-6))\V_hat'*Ry_hat_inv*Y);
        if (B_now == B_hat) | (d==threshold)
            break;
        else
            B_hat = B_now;
        end
    end
    extractedWatermark = B_hat;
end
