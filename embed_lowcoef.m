function [Y,Ry_inv,Ry_hat_inv,s] = embed_lowcoef(xm,B,alpha,noise,sourceNum,len)
% Embed B into host.
% parameters: input: host:xm; watermark:B; amplitude:alpha; 
%                    noise power:noise; number of carriers:sourceNum;
%                    length of embedded coefficients(from bottom to top):len
% output: embedded matrix:Y; inverse of Y's autocorrelation matrix:Y; inverse of sample average autocorrelation matrix:Ry_hat_inv
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
    [a,b] = size(xm);
    
   if sourceNum == 8 && len == 8
        %U is generated randomly
        U = [-0.394428407449048,-0.0421792326291786,-0.0340380082682522,-0.146850031139384,-0.0916921833470477,0.153883110399089,-0.144689592162555,-0.538439624162039;-0.488091245481685,-0.592892092324800,0.290861073265807,-0.263071186016667,0.816526154134713,-0.361152536933426,-0.308523991755411,0.140529959343919;-0.583631605586807,0.289509741777581,-0.205135226437268,0.519159947442117,-0.457331690268016,-0.602635798275545,-0.208102081252585,0.0989961687961332;0.225532270284864,0.268400475509108,-0.167448292833135,-0.136402501516466,-0.00300567845524254,-0.102307857312898,-0.321516020066217,0.654586333402916;0.0373753948395372,-0.387098572362618,0.400158425804600,0.253333525446319,0.197189401480857,0.280732410739488,-0.247112043670124,-0.144781470493184;-0.0603740517708346,0.415966549796958,0.182061342858689,0.383321441028059,0.0345062746024347,-0.288797823599255,-0.0441149113472229,0.246977180115821;0.454140445282694,0.228062702693872,-0.803481985050361,0.477470691277013,0.274868934861799,-0.370924105361524,-0.692738203325451,0.156544747570853;0.0584552922622377,-0.340527681284066,0.0725634863219728,0.426619810485469,-0.00943396551145083,0.415273968695094,-0.440797158763178,-0.381555936471573];
        s = U;
        U = U*alpha;
  
   else
        %make sure the function works if not 8x8 lattice
        U = zeros(len,sourceNum);
        for i = 1:sourceNum
            U(:,i) = randn(len,1);
            U(:,i) = U(:,i)/norm(U(:,i));
        end
        s = U;
        U=U*alpha;  
    end
    
    %calculate Rx
    Rx=zeros(len,len);
    for i=1:b
        tmp=xm(a-len+1:a,i);
        Rx=Rx+tmp*tmp';
    end
    Rx=Rx/b;
    
    %generate N
    if noise == 0
        N = zeros(len,b);
    else
        N = zeros(len,b);
        for i = 1:b
           N(:,i) = randn(len,1)*noise; 
        end
    end
    
    %embed
    Z = xm(a-len+1:a,:) + N;
    Y = U*B + Z;
    
    %calculate Ry_inv
    Ry = Rx+eye(len,len)*(noise^2)+U*U';
    Ry_inv = inv(Ry);
    
    %calculate Ry_hat_inv
    Ry_hat = zeros(len,len);
    for i = 1:b
        tmp = Y(:,i);
        Ry_hat = Ry_hat + tmp*tmp';
    end
    Ry_hat = Ry_hat / b;
    Ry_hat_inv = inv(Ry_hat);
end