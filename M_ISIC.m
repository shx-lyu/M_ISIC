function [BERs,iterate_times] = M_ISIC(Y,B,Ry_hat_inv,P,threshold,sourceNum,host_height,host_width,Bc_hatsic)
% Do P times M_ISIC,compute P bers
% Parameters: input: embedded matrix:Y; prewhitening matrix:Ry_hat_inv;
%                    origin watermark:B; reinitialization times:P;
%                    threshold of iterating times:threshold;
%                    number of carriers:sourceNum
%                    pre-initialization:Bc_hatsic
% Output: 1xP vector of BER:BERs; M_ISIC's actual iterating times:iterate_times
k = 1;
BERs = zeros(1,P);
iterate_times = zeros(1,P);
host_hig_wid = host_height*host_width;
Y = sqrtm(Ry_hat_inv)*Y;                                    % prewhitening,Y=Ry_inv_sqr*Y
V_here=Y*Bc_hatsic{1,1}'/(Bc_hatsic{1,1}*Bc_hatsic{1,1}');  
[M,N]=size(V_here);
while(k<=P)
    [extractedWatermark,iterate_times(k)] = SIC_extract(Y,threshold,sourceNum,host_hig_wid,Bc_hatsic{1,k},M,N);
    [BERs(k),~] = BER(B,extractedWatermark,sourceNum,host_height,host_width);
    k = k+1;
end
end

function [extractedWatermark,d] = SIC_extract(Y,threshold,sourceNum,host_hig_wid,B_hat,M,N)
d = 0;
cons=[-1 1];
one_1=cons(1)*ones(1,host_hig_wid); 
one_2=cons(2)*ones(1,host_hig_wid);
B_now = zeros(sourceNum,host_hig_wid);

while(1)
    d = d+1;
    V_hat = Y*B_hat'/(B_hat*B_hat');
    
% ------------- sorting ----------------- %
    len_H=diag(V_hat'*V_hat);
    [~,IND]=sort(len_H);
    Vsorted=V_hat(:,IND);
% ------------- sorting ----------------- %

    [Q,R2] = qr(Vsorted);
    Y_q=Q(1:M,1:N)'*Y;
    R=R2(1:N,1:N);
    
    for  j=N:-1:1                   
        tmp_y_row=Y_q(j,:)/R(j,j);  
        cons_2=abs([one_1-tmp_y_row ; one_2-tmp_y_row]);
        B_now(j,:)=sign(cons_2(1,:)-cons_2(2,:));
        if j>=2                     
            Y_q=Y_q-R(:,j)*B_now(j,:);
        end
    end
    
    B_now(IND,:)=B_now(1:N,:);
    if (norm(B_now-B_hat,'fro')<1e-5) | (d==threshold)
        break;
    else
        B_hat = B_now;
    end
end
extractedWatermark = B_hat;
end