function v = zigzag(a)
% Do a zigzag scan on an 8x8 block,only 8x8 will work
% parameters: input: 8x8 matrix:a
% output: 1x64 vector v
% author: Hao Cheng, haoc678@gmail.com, Jinan University
%         Fan Yang, fanyang1124@163.com, Jinan University
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
