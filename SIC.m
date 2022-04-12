% lattice-based alternating minimization+Babai(SIC) method
% ref: Shanxiang Lyu, Zheng Wang, Zhen Gao, Hongliang He, Lajos Hanzo, "Lattice-Based mmWave Hybrid Beamforming", IEEE Transactions on Communications,2021
% Author  : Shanxiang Lyu (shanxianglyu@gmail.com, https://sites.google.com/view/shanx)
% Date    : 2020-June
function [x_hat]=SIC(yy,H,Q,R2,cons)

    [M,N]=size(H);
%     [Q,R2]=qr(H); 

    if M>N
        Q1=Q(1:M,1:N);
        y=Q1'*yy;
    else
        y=Q'*yy;
    end
    R=R2(1:N,1:N); 

    x_hat=zeros(N,1);
    for i=N:-1:1
        x_hat(i)=ClimbOne(y(i)/R(i,i),cons);
        if i>=2
        y=y-R(:,i)*x_hat(i);
        end
    end
end

function s_out=ClimbOne(s_in,cons)

len=size(cons,2);

s_com2=ones(1,len)*s_in;

s_com=abs(cons-s_com2);

if s_com(1)<s_com(2)
    s_out = -1;
else
    s_out = 1;
end

end