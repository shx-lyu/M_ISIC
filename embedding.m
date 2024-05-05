function [Y,Ry_inv,Ry_hat_inv,s,U] = embedding(xm,B,alpha,noise,U,len)
% Embed B into host.
% Parameters: input: host:xm; watermark:B; amplitude:alpha; 
%                    noise power:noise; number of carriers:sourceNum;
%                    length of embedded coefficients(from bottom to top):len
% Output: embedded matrix:Y; inverse of Y's autocorrelation matrix:Y; inverse of sample average autocorrelation matrix:Ry_hat_inv
    [a,b] = size(xm);
    s = U;
    U = U*alpha;
    
    %calculate Rx
    Rx=zeros(len,len);
    for i=1:b
        tmp=xm(a-len+1:a,i);
        Rx=Rx+tmp*tmp';
    end
    Rx=Rx/b;
    
    % generate N
    if noise == 0
        N = zeros(len,b);
    else
        N = zeros(len,b);
        for i = 1:b
           N(:,i) = randn(len,1)*noise; 
        end
    end
    
    % embed
    Z = xm(a-len+1:a,:) + N;
    Y = U*B + Z;
    
    % calculate Ry_inv
    Ry = Rx+eye(len,len)*(noise^2)+U*U';
    Ry_inv = inv(Ry);
    
    % calculate Ry_hat_inv
    Ry_hat = zeros(len,len);
    for i = 1:b
        tmp = Y(:,i);
        Ry_hat = Ry_hat + tmp*tmp';
    end
    Ry_hat = Ry_hat / b;
    Ry_hat_inv = inv(Ry_hat);
end