% Non-Blind Extraction: Comparing the performance of Sphere Decoding and other algrithms.
clc;
clear;
tic;
% --------------------------------------- Parameters --------------------------------------- %
Distortion = [24,26,28,30,32,34,36,38];  % Distortion per message
alphas = sqrt(10.^(Distortion/10));      % A_k
noise = sqrt(3);                         % Noise
sourceNum = 8;                           % Number of carriers
len = 8;                                 % Dimesion of the lattice
N = 20;                                  % Simulation times
% ------------------------------------------------------------------------------------------ %
host = 'pic/b.jpg';
[xm,Rx,dc,host_height,host_width] = DctHost(host);
host_height=host_height/8;
host_width=host_width/8;
w_len=host_height*host_width;

ber_mmse = zeros(1, length(alphas));
ber_zf = zeros(1, length(alphas));
ber_sic = zeros(1, length(alphas));
ber_sd = zeros(1, length(alphas));
ber_ideal = zeros(1, length(alphas));

Mat_Gen = Matrix_Generation;
Uc = Mat_Gen.Gen_U(N,sourceNum,len);
Bc = Mat_Gen.Gen_B(N,sourceNum,w_len);
MMSE = MMSE_filters;
ALA = Approx_Lat_Algo;
for i = 1:length(alphas)
    disp('i= '+string(i));
    n_mmse = 0;
    n_zf = 0;
    n_sic = 0;
    n_sd = 0;
    n_ideal = 0;
    for j = 1:N
        disp('Distortion=   '+string(Distortion(i))+'   image=   '+ host+'   k=   '+ string(j));
        [Y,Ry_inv,Ry_hat_inv,s,U] = embedding(xm,Bc{1,j},alphas(i),noise,Uc{1,j},len);
      
        % SMI-MMSE
        n_mmse = n_mmse + MMSE.SMI_MMSE(s,Ry_hat_inv,Y,Bc{1,j},sourceNum,host_height,host_width);
        disp('   SMI-MMSE   Done!');
        
        % zero-forcing
        n_zf = n_zf + ALA.Zero_Forcing(Y, U, Bc{1,j}, Ry_hat_inv);
        disp('   Zero Forcing   Done!');
        
        % sic
        n_sic = n_sic + ALA.SIC(Y, U, Bc{1,j}, Ry_hat_inv);
        disp('   SIC   Done!');
        
        % SD
        n_sd = n_sd + sph_dec(Y, U, Bc{1,j}, Ry_hat_inv);
        disp('   Sphere Decoding   Done!');
        
        % Ideal-MMSE
        n_ideal = n_ideal + MMSE.Ideal_MMSE(s,Ry_inv,Y,Bc{1,j},sourceNum,host_height,host_width);
        disp('   Ideal-MMSE   Done!');
    end
    ber_mmse(i) = n_mmse / N;
    ber_zf(i) = n_zf / N;
    ber_sic(i) = n_sic / N;
    ber_sd(i) = n_sd / N;
    ber_ideal(i) = n_ideal / N;
end

x = Distortion;
figure
axes('yscale', 'log')
hold on
xlabel('Distortion in dB (per-message)');
hold on
ylabel('Average BER');
hold on
semilogy(x,ber_zf,'r-^',x,ber_sic,'b-',x,ber_mmse,'k--', x,ber_sd,'m-d',x,ber_ideal,'g-^');

% semilogy(x,BER_matrix(1,:),'r-^',x,BER_matrix(2,:),'k-',x,BER_matrix(3,:),'k--',x,BER_matrix(4,:),'r-d',x,BER_matrix(5,:),'k.-',x,BER_matrix(6,:),'g-');
grid on
legend({'ZF','SIC','SMI-MMSE','SD','Ideal-MMSE'},'Location','southwest');
toc;







