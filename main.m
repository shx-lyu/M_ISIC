% main function:Comparing the performance of M_ISIC and other extracting algrithms.
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University

clc;
clear;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Distortion = [24,26,28,30,32,34,36,38]; %Distortion per message

alphas = sqrt(10.^(Distortion/10)); %A_k

noise = sqrt(3);    %Noise

sourceNum = 8;      %Number of carriers  

len = 8;            %Dimesion of the lattice

threshold = 300;    %Threshold of iterating times of M_IGLS and M_ISIC,make sure the algorithm stops.

P = 20;             %Reinitialization times of M_IGLS and M_ISIC

Q = 1;              %Image numbers

N = 20;              %Simulation times

range = 0.5;        %interval of quantizing to 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BER_matrix = zeros(6,length(alphas));

for i = 1:length(alphas)
    
    %init
    M_IGLS_sum = 0;
    SMI_MMSE_sum = 0;
    Ideal_MMSE_sum = 0;
    JADE_sum = 0;
    MF_sum = 0;
    SIC_sum = 0;
    
    for j = 1:Q
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If you want to run this code over a set of images,
        % please replace the code below with a loop
        
        % host = '52.pgm';  %A low information entropy image
        host = 'gray.png';  %A high information entropy image
        [xm,Rx,dc,host_height,host_width] = DctHost(host);
        
        % After the following 2 lines of code are executed,
        % 'host_height' and 'host_width' are no longer the actual size of the host
        % image. 'host_height*host_height' represents the sum of the blocks.
        host_height=host_height/8;
        host_width=host_width/8;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('Distortion=   '+string(Distortion(i))+'   image=   '+string(j));
        
        %Pre-initialize Bc_hat£¬make sure M_IGLS and M_ISIC have same initialization
        Bc_hat=cell(1,P);
        Bc_hatsic=cell(1,P);
        for p=1:P
            Bc_hat{1,p}=randn(sourceNum,host_height*host_width);
            Bc_hat{1,p}=sign(Bc_hat{1,p});
            Bc_hatsic{1,p}= Bc_hat{1,p};
        end
        
        %M-IGLS
        for k = 1:N
            B = watermarks2B(sourceNum,host_height,host_width);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            [IGLS_tmp_bers,d_IGLS] = M_IGLS(Y,Ry_hat_inv,B,P,threshold,sourceNum,host_height,host_width,Bc_hat);
            M_IGLS_sum = M_IGLS_sum+min(IGLS_tmp_bers);
        end
        disp('   M-IGLS     Done!');
        
        %SMI-MMSE
        for k = 1:N
            B = watermarks2B(sourceNum,host_height,host_width);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            SMI_MMSE_sum = SMI_MMSE_sum+SMI_MMSE(s,Ry_hat_inv,Y,B,sourceNum,host_height,host_width);
        end
        disp('   SMI-MMSE   Done!');
        
        %Ideal-MMSE
        for k = 1:N
            B = watermarks2B(sourceNum,host_height,host_width);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            Ideal_MMSE_sum = Ideal_MMSE_sum+Ideal_MMSE(s,Ry_inv,Y,B,sourceNum,host_height,host_width);
        end
        disp('   Ideal-MMSE Done!');
        
        %JADE
        for k = 1:N
            B = watermarks2B(sourceNum,host_height,host_width);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            JADE_sum = JADE_sum+JADE(Y,sourceNum,B,host_height,host_width);
        end
        disp('   JADE       Done!');
        
        %MF
        for k = 1:N
            B = watermarks2B(sourceNum,host_height,host_width);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            [MF_ber,~] = BER(B,sign(s'*Y),sourceNum,host_height,host_width);
            MF_sum = MF_sum+MF_ber;
        end
        disp('   MF         Done!');
        
        %SIC
        for k = 1:N
            B = watermarks2B(sourceNum,host_height,host_width);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            [SIC_tmp_bers,d_SIC] = Lattice_Based_extract(Y,B,Ry_hat_inv,P,threshold,sourceNum,host_height,host_width,Bc_hatsic);
            SIC_sum = SIC_sum+min(SIC_tmp_bers);
        end
        disp('   SIC        Done!');
 
    end
    
    %calculate average BER
    BER_matrix(1,i) = M_IGLS_sum/(Q*N);
    BER_matrix(2,i) = SMI_MMSE_sum/(Q*N);
    BER_matrix(3,i) = Ideal_MMSE_sum/(Q*N);
    BER_matrix(4,i) = JADE_sum/(Q*N);
    BER_matrix(5,i) = MF_sum/(Q*N);
    BER_matrix(6,i) = SIC_sum/(Q*N);
end
disp('Done!');

%plot
x = Distortion;
figure
axes('yscale', 'log')
hold on
xlabel('Distortion in dB (per-message)');
hold on
ylabel('Average BER');
hold on
semilogy(x,BER_matrix(1,:),'r-^',x,BER_matrix(2,:),'k-',x,BER_matrix(3,:),'k--',x,BER_matrix(4,:),'r-d',x,BER_matrix(5,:),'k.-',x,BER_matrix(6,:),'g-');
grid on
legend({'M-IGLS','SMI-MMSE','Ideal-MMSE','JADE','MF','SIC'},'Location','southwest');

toc;