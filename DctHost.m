function [xm,Rx,dc,host_height,host_width] = DctHost(host)
% Do DCT transform on the host iamge
% Parameters: input:host iamge path:host
% Output: DCT coefficients:xm; Autocorrelation matrix of host image:Rx;
%         DC component:dc; size of host image:host_height,host_width
    host = imread(host);
    host = im2double(host)*255;
    [host_height,host_width] = size(host);
    
    % Do DCT transform on every 8x8 block 
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
    
    % calculate autocorrelation matrix of host image
    Rx = zeros(63,63);
    for i = 1:(host_height/8)*(host_width/8)
        tmp = xm(:,i);
        Rx = Rx + tmp*tmp';
    end
    Rx = Rx/((host_height/8)*(host_width/8));
end

function v = zigzag(a)
% Do a zigzag scan on an 8x8 block,only 8x8 will work
% Parameters: input: 8x8 matrix:a
% Output: 1x64 vector v
   [n,m] = size(a);
   v = zeros([1 64]);
   if(n~=8||m~=8)
      disp('Not a 8x8 matrix!'); 
      return
   end
   i = 1;   %row
   j = 1;   %column
   k = 1;
   l = 1;   %iteration times
   times = [1 2 3 4 5 6 7 6 5 4 3 2 1];
   while(i~=9 && j~=9)
        v(k) = a(i,j);
        k = k+1;

        if(i==8 && j==7)
            break
        elseif(i==1)
            j = j+1;
            for time = 1:times(l)
                v(k) = a(i,j);
                i = i+1;
                j = j-1;
                k = k+1;
            end
            l = l+1;

        elseif(j==1 && i~=8)
            i = i+1;
            for time = 1:times(l)
                v(k) = a(i,j);
                i = i-1;
                j = j+1;
                k = k+1;
            end
            l = l+1;

        elseif(i==8)
            j = j+1;
            for time = 1:times(l)
                v(k) = a(i,j);
                i = i-1;
                j = j+1;
                k = k+1;
            end
            l = l+1;

        elseif(j==8)
            i = i+1;
            for time = 1:times(l)
                v(k) = a(i,j);
                i = i+1;
                j = j-1;
                k = k+1;
            end
            l = l+1;
        end
   end
   v(63) = a(8,7);
   v(64) = a(8,8);
end
