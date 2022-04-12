function [BERs,iterate_times] = Lattice_Based_extract(Y,B,Ry_hat_inv,P,threshold,sourceNum,host_height,host_width,Bc_hatsic)
% Do P times M_ISIC,compute P bers
% parameters: input: embedded matrix:Y; prewhitening matrix:Ry_hat_inv;
%                    origin watermark:B; reinitialization times:P;
%                    threshold of iterating times:threshold;
%                    number of carriers:sourceNum
%                    pre-initialization:Bc_hatsic
% output: 1xP vector of BER:BERs; M_ISIC's actual iterating times:iterate_times
% author: Shanxiang Lyu (shanxianglyu@gmail.com, https://sites.google.com/view/shanx)
%         Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
k = 1;
BERs = zeros(1,P);
iterate_times = zeros(1,P);
host_hig_wid = host_height*host_width;
Y = sqrtm(Ry_hat_inv)*Y;                                    % prewhitening,Y=Ry_inv_sqr*Y
V_here=Y*Bc_hatsic{1,1}'/(Bc_hatsic{1,1}*Bc_hatsic{1,1}');  % gain size of V_here(V_here is temporary variable)
[M,N]=size(V_here);
while(k<=P)
    [extractedWatermark,iterate_times(k)] = SIC_extract(Y,threshold,sourceNum,host_hig_wid,Bc_hatsic{1,k},M,N);
    [BERs(k),~] = BER(B,extractedWatermark,sourceNum,host_height,host_width);
    k = k+1;
end
end

function [extractedWatermark,d] = SIC_extract(Y,threshold,sourceNum,host_hig_wid,B_hat,M,N)
% Core code of M_ISIC
d = 0;
cons=[-1 1];
one_1=cons(1)*ones(sourceNum,host_hig_wid); 
one_2=cons(2)*ones(sourceNum,host_hig_wid);

while(1)
    d = d+1;
    
    V_hat = (Y*B_hat')/(B_hat*B_hat');
    len_H = diag(V_hat'*V_hat);
    [~,IND] = sort(len_H);
    Vsorted = V_hat(:,IND);    
    
    [Q,R2] = qr(Vsorted);       % V_hat = (Y*B_hat')/(B_hat*B_hat')
    
%     if M>N
        Y_q = Q(1:M,1:N)'*Y;                      % run this
%     else
%         Y_q = Q'*Y;                 
%     end
%    R=R2(1:N,1:N);  
    digR=repmat(diag(R2),1,host_hig_wid);
    Y_q=Y_q./digR;
    B_now=sign(abs(one_1-Y_q)-abs(one_2-Y_q));
    
    if (B_now == B_hat) | (d==threshold)
        break;
    else
        B_hat = B_now;
    end
end
extractedWatermark = B_hat;
end