function r = SDCVP(z, H, R)
% return the x coefficient, r; z = Q'*y (对应论文y=Q1^{*}x)
radius = realmax;
n = size(H,2);

% add examine this variable before make it global
global SPHDEC_RADIUS;
global RETVAL;
global x;
global NUMX;

SPHDEC_RADIUS = radius;  % 搜索半径会动态调整, 初始为无穷大
RETVAL = zeros(n, 1);    % n 表示 H 的列数
x = zeros(n, 1);
NUMX  = 0;

sphdec_core(z, R, n, 0);
if NUMX > 0
    r = RETVAL;          % retrival is set to r only this time
else
    r = zeros(n, 1);     % no vector inside, return zeros
end
clear SPHDEC_RADIUS RETVAL SYMBSETSIZE NUMX;

% z = Q'*y
function sphdec_core(y_bar, R, layer, dist)
% dist and d are both current radius (使用全局变量时也要用global定义)
global SPHDEC_RADIUS;
global RETVAL;
global x;
global NUMX;

n=size(R,2);
if layer==n                                           % 如果是最后一层
    y_bpr=y_bar;                                      % the last y is the first y
else
    y_bpr=y_bar - R(:,layer+1:end)*x(layer+1:end);    % Q1^{T}x-Rs, minius coef in prev rounds
end

c = sign(y_bpr(layer)/R(layer,layer));
C = [c, -c];
if (layer == 1)
    for i = 1:2
        x(layer) = C(i);
        d = abs(y_bar(layer) - R(layer,layer:end)*x(layer:end))^2 + dist;
        if  (d <= SPHDEC_RADIUS)
            RETVAL        =  x;     % 保存x
            SPHDEC_RADIUS =  d;     % 如果找到一个可行解, 则缩小半径(初始时半径为无穷)
            NUMX   =  NUMX + 1;     % 可行解+1
        end
    end
else
    for i = 1:2
        x(layer) = C(i);
        d = abs(y_bar(layer) - R(layer,layer:end)*x(layer:end))^2 + dist;
        if  (d <= SPHDEC_RADIUS)
            sphdec_core(y_bar, R,  layer-1, d);
        end
    end  
end

% for i = 1:2
%     x(layer) = C(i);
%     d = abs(y_bar(layer) - R(layer,layer:end)*x(layer:end))^2 + dist;
%     if (layer == 1) && (d <= SPHDEC_RADIUS)
%             RETVAL        =  x;     % 保存x
%             SPHDEC_RADIUS =  d;     % 如果找到一个可行解, 则缩小半径(初始时半径为无穷)
%             NUMX   =  NUMX + 1;     % 可行解+1
%     end
%     if (layer > 1) && (d <= SPHDEC_RADIUS)
%             sphdec_core(y_bar, R,  layer-1, d);
%     end
% end
        

