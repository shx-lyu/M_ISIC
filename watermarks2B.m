function B = watermarks2B(sourceNum,host_height,host_width)
% Generate a random column iid. watermark NxM matrix B\in \{\pm 1\},N=sourceNum,M=host_height*host_width
% Parameters: input: number of carriers:sourceNum;
% Output: a random watermark matrix:B
% Author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
     B = zeros(sourceNum,host_height*host_width);
     for i = 1:host_height*host_width
         B(:,i) = sign(randn(sourceNum,1));
     end
end