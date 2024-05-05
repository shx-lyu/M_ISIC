function [ber,correctWatermark] = BER(B,extractedWatermark,sourceNum,host_height,host_width)
% Compute BER of extractedWatermark using greedy method, and then correct its order and sign(since blind extract algorithm causes sign-ambiguity problem)
% Parameters: input: orgin watermark:B; extracted watermark:extractedWatermark
%                    number of carriers:sourceNum
% Output: bit error rate of extractedWatermark:ber; 
%         watermark which order and sign are corrected:correctWatermark
    correctWatermark = zeros(sourceNum,host_height*host_width);
    for i = 1:sourceNum
        tmp = extractedWatermark(i,:);
        minRate = 2;
        flag = 1;
        position = 1;
        for j = 1:sourceNum
           if correctWatermark(j,1) ~= 0
               continue
           end
           tmp1 = (tmp == B(j,:));
           tmp2 = (tmp == (-1)*B(j,:));
           rate1 = length(find(tmp1==0))/(host_height*host_width);
           rate2 = length(find(tmp2==0))/(host_height*host_width);
           if rate1<minRate
               minRate = rate1;
               position = j;
               flag = 1;
           end
           if rate2<minRate
               minRate = rate2;
               position = j;
               flag = -1;
           end
        end
        correctWatermark(position,:) = flag*tmp;
    end
    tmp = (correctWatermark == B);
    ber = length(find(tmp==0))/(sourceNum*(host_height*host_width));
end