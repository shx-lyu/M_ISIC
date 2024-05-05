% Blind Extraction: Comparing the performance of M_ISIC and other algrithms.
clc;
clear;
tic;
% ------------------------------------- Parameters ------------------------------------- %
Distortion = [24,26,28,30,32,34,36,38] ; % Distortion per message
alphas = sqrt(10.^(Distortion/10));      % A_k
noise = sqrt(3);     % Noise
sourceNum = 8;       % Number of carriers  
len = 8;             % Dimesion of the lattice
threshold = 50;      % Threshold of iterating times
P = 20;              % Reinitialization times of M_IGLS and M_ISIC
N = 20;              % Simulation times
% ------------------------------------------------------------------------------------- %
host = 'pic/gray.png';  
[xm,Rx,dc,host_height,host_width] = DctHost(host);
host_height=host_height/8;
host_width=host_width/8;
w_len=host_height*host_width;

Mat_Gen = Matrix_Generation;
Uc = Mat_Gen.Gen_U(N,sourceNum,len);
Bc_hat = Mat_Gen.Gen_B(P,sourceNum,w_len);
Bc_watermarks2B = Mat_Gen.Gen_B(N,sourceNum,w_len);

BER_matrix = zeros(5,length(alphas));
MMSE = MMSE_filters;
for i = 1:length(alphas)
    
    M_IGLS_sum = 0;
    SMI_MMSE_sum = 0;
    Ideal_MMSE_sum = 0;
    JADE_sum = 0;
    SIC_sum = 0;
    for k = 1:N
        disp('Distortion=   '+string(Distortion(i))+'   image=   '+ host+'   k=   '+ string(k));
        [Y,Ry_inv,Ry_hat_inv,s,~] = embedding(xm,Bc_watermarks2B{1,k},alphas(i),noise,Uc{1,k},len);
            
        % M-IGLS
        [IGLS_tmp_bers,d_IGLS] = M_IGLS(Y,Ry_hat_inv,Bc_watermarks2B{1,k},P,threshold,sourceNum,host_height,host_width,Bc_hat);
        M_IGLS_sum = M_IGLS_sum + min(IGLS_tmp_bers);
        disp('   M-IGLS     Done!');
        
        % SMI-MMSE
        SMI_MMSE_sum = SMI_MMSE_sum + MMSE.SMI_MMSE(s,Ry_hat_inv,Y,Bc_watermarks2B{1,k},sourceNum,host_height,host_width);
        disp('   SMI-MMSE   Done!');
        
        % Ideal-MMSE
        Ideal_MMSE_sum = Ideal_MMSE_sum + MMSE.Ideal_MMSE(s,Ry_inv,Y,Bc_watermarks2B{1,k},sourceNum,host_height,host_width);
        disp('   Ideal-MMSE Done!');
        
        % JADE
        JADE_sum = JADE_sum+JADE(Y,sourceNum,Bc_watermarks2B{1,k},host_height,host_width);
        disp('   JADE       Done!');
            
        % M-ISIC
        [SIC_tmp_bers,d_SIC] = M_ISIC(Y,Bc_watermarks2B{1,k},Ry_hat_inv,P,threshold,sourceNum,host_height,host_width,Bc_hat);
        SIC_sum = SIC_sum+min(SIC_tmp_bers);
        disp('   M-ISIC    Done!');
    end
 
    % average BER
    BER_matrix(1,i) = M_IGLS_sum/N;
    BER_matrix(2,i) = SMI_MMSE_sum/N;
    BER_matrix(3,i) = Ideal_MMSE_sum/N;
    BER_matrix(4,i) = JADE_sum/N;
    BER_matrix(5,i) = SIC_sum/N;
end
disp('Done!');

% plot
x = Distortion;
figure
axes('yscale', 'log')
hold on
xlabel('Distortion in dB (per-message)');
hold on
ylabel('Average BER');
hold on
semilogy(x,BER_matrix(1,:),'r-^',x,BER_matrix(2,:),'k-',x,BER_matrix(3,:),'k--',x,BER_matrix(4,:),'r-d',x,BER_matrix(5,:),'g-');
grid on
legend({'M-IGLS','SMI-MMSE','Ideal-MMSE','JADE','M-ISIC'},'Location','southwest');

toc;