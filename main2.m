clc;
clear;
tic;

% host = '52.pgm';   %A low information entropy image
host = 'gray.png';   %A high information entropy image
[xm,Rx,dc,host_height,host_width] = DctHost(host);
host_height=host_height/8;
host_width=host_width/8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Distortion = [24,26,28,30,32,34,36,38];
alphas = sqrt(10.^(Distortion/10));

noise = sqrt(3);    %Noise

sourceNum = 4;      %Number of carriers  

len = 10;           %Dimesion of the lattice

threshold = 300;    %Threshold of iterating times of M_IGLS and M_ISIC

guess = 5;          %Number of guesses of carriers

P = 20;             %Reinitialization times of M_IGLS and M_ISIC

Q = 1;              %Image numbers

N = 1;              %Simulation times

range = 0.5;        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%存放每一种算法的ber的平均值

BER_matrix = zeros(4,length(Distortion));

%对每一个alpha做Q次嵌入，提取后取平均
%对于M-IGLS和SIC，做P次二次迭代然后找到最小的BER
%对于SMI_MMSE，Ideal_MMSE，JADE，FastICA，每一张图片提取后累加，最后算平均值
for i = 1:length(Distortion)
    
    %L = 5 6 7 8
    for j = 1:length(guess)
        
        Ideal_MMSE_sum = 0;
        SMI_MMSE_sum = 0;
        SIC_sum = 0;
        IGLS_sum = 0;
        
        disp('L   =   '+string(guess(j))+'   Distortion   =   '+string(Distortion(i)));
        
        %预先对Bc_hat初始化，保证gls和sic初始化相同
        Bc_hat=cell(1,P);
        Bc_hatsic=cell(1,P);
        for p=1:P
            Bc_hat{1,p}=randn(guess(j),host_height*host_width);
            Bc_hat{1,p}=sign(Bc_hat{1,p});
        end
        
        for p = 1:P
           Bc_hatsic{1,p}=randn(guess(j),host_height*host_width);
           Bc_hatsic{1,p}=sign(Bc_hatsic{1,p});
        end
        
        %Ideal_MMSE
        for k = 1:N
            B = watermarks2B(sourceNum,64,64);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            Ideal_MMSE_sum = Ideal_MMSE_sum+Ideal_MMSE(s,Ry_inv,Y,B,sourceNum,host_height,host_width);
        end
        disp('   Ideal-MMSE   Done!');
        BER_matrix(1,i) = Ideal_MMSE_sum/(Q*N); 
        
       %SMI-MMSE
        for k = 1:N
            B = watermarks2B(sourceNum,64,64);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            SMI_MMSE_sum = SMI_MMSE_sum+SMI_MMSE(s,Ry_hat_inv,Y,B,sourceNum,host_height,host_width);
        end
        disp('   SMI-MMSE   Done!');
        BER_matrix(2,i) = SMI_MMSE_sum/(Q*N);
        
        %SIC
        for k = 1:N
            B = watermarks2B(sourceNum,64,64);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            [SIC_tmp_bers,d_SIC] = Lattice_Based_extract_new(Y,B,eye(len,len),P,threshold,sourceNum,host_height,host_width,Bc_hatsic,guess(j),range);
%             [SIC_tmp_bers,d_SIC] = Lattice_Based_extract_new(Y,B,Ry_hat_inv,P,threshold,sourceNum,host_height,host_width,Bc_hatsic,guess(j),range);
            SIC_sum = SIC_sum+min(SIC_tmp_bers);
        end
        disp('   SIC        Done!');
        BER_matrix(3,i) = SIC_sum/N; 
        
        %IGLS
        for k = 1:N
            B = watermarks2B(sourceNum,64,64);
            [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alphas(i),noise,sourceNum,len);
            [IGLS_tmp_bers,d_IGLS] = M_IGLS(Y,eye(len,len),B,P,threshold,sourceNum,host_height,host_width,Bc_hat,guess(j));
            IGLS_sum = IGLS_sum+min(IGLS_tmp_bers);
        end
        disp('   IGLS        Done!');
        BER_matrix(4,i) = IGLS_sum/N;
    end
end
disp('Done!');

%画图
x = Distortion;
figure
axes('yscale', 'log')
hold on
xlabel('Distortion in dB (per message)');
hold on
ylabel('Average BER');
hold on
% ylim([10e-4,10e-3]);
semilogy(x,BER_matrix(1,:),'b--',x,BER_matrix(2,:),'b.-',x,BER_matrix(3,:),'r^-',x,BER_matrix(4,:),'r-o');
grid on
legend({'Ideal-MMSE','SMI-MMSE','M-ISIC','M-IGLS'},'Location','southwest');
toc;