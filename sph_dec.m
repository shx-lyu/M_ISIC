function [error_rate] = sph_dec(Y, V, B, Ry_hat_inv)
Y = sqrtm(Ry_hat_inv) * Y;
[Q, R] = qr(sqrtm(Ry_hat_inv) * V,0);

Y_bar = Q'*Y;
[m, n] = size(B);
B_hat = zeros(m, n);

for i = 1:n
    y_bar = Y_bar(:, i);
    B_hat(:, i) = SDCVP(y_bar, V, R);
end
ber_mat = (B_hat == B);
error_rate = length(find(ber_mat==0))/(size(B,1)*size(B,2));
end

