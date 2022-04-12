function [xm,Rx,dc,host_height,host_width] = DctHost(host)
% Do DCT transform on the host iamge
% parameters: input:host iamge path:host
% output: DCT coefficients:xm; Autocorrelation matrix of host image:Rx;
%         DC component:dc; size of host image:host_height,host_width
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
    host = imread(host);
    host = im2double(host)*255;
    [host_height,host_width] = size(host);
    
    %Do DCT transform on every 8x8 block 
    T = dctmtx(8);
    dct = @(block_struct) T*block_struct.data*T';
    host_dct = blockproc(host,[8 8],dct);
    zigzagscan = @(block_struct) zigzag(block_struct.data);
    x = blockproc(host_dct,[8 8],zigzagscan);
    
    xm = zeros(64,(host_height/8)*(host_width/8));
    for i = 1:host_height/8
        for j = 1:host_width/8
            xm(:,(i-1)*(host_height/8)+j) = x(i,(j-1)*64+1:j*64);
        end
    end
    
    dc = xm(1,:);
    [~,w] = size(xm);
    tmp = zeros(63,w)+xm(2:64,:);
    xm = tmp;
    
    %calculate autocorrelation matrix of host image
    Rx = zeros(63,63);
    for i = 1:(host_height/8)*(host_width/8)
        tmp = xm(:,i);
        Rx = Rx + tmp*tmp';
    end
    Rx = Rx/((host_height/8)*(host_width/8));
end