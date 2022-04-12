% 计算多维噪声的最大值
% 假设二维噪声，iid，A_1^2+A_2^2=1
clc;
clear;
A1_square=0.01:0.01:1;
A2_square=1-A1_square;
integral1=normcdf(sqrt(A1_square))-normcdf(-sqrt(A1_square));
integral2=normcdf(sqrt(A2_square))-normcdf(-sqrt(A2_square));
d=integral1.*integral2;

figure(1);
grid on;
subplot(1,2,1);
plot(A1_square,A2_square);

grid on;
subplot(1,2,2);
plot3(A1_square,A2_square,d);
grid on;
