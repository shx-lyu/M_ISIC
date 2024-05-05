% Low complexity lattice decoding algorithms: SIC and Zero Forcing
function ALA = Approx_Lat_Algo
ALA.SIC = @SIC;
ALA.Zero_Forcing = @Zero_Forcing;
end

function [error_rate] = SIC(Y, V, B, Ry_hat_inv)
Y = sqrtm(Ry_hat_inv)*Y;
cons = [-1 1];
[m, n] = size(B);
B_hat = zeros(m, n);
one_1 = cons(1)*ones(1, n);
one_2 = cons(2)*ones(1, n);

% QR decomposition
[Q, R2] = qr(sqrtm(Ry_hat_inv)*V);
[M, N] = size(V);

if M > N
    Y_q = Q(1:M, 1:N)'*Y;
else
    Y_q = Q'*Y;
end
R = R2(1:N,1:N);

for j = N:-1:1 
    tmp = Y_q(j,:)/R(j,j);
    cons_2 = abs([one_1 - tmp; one_2 - tmp]);
    B_hat(j,:) = sign(cons_2(1,:)-cons_2(2,:));
    if j >= 2
        Y_q = Y_q-R(:,j)*B_hat(j,:);
    end
end

ber_mat = (B_hat == B);
error_rate = length(find(ber_mat==0))/(size(B,1)*size(B,2));
end

function [error_rate] = Zero_Forcing(Y, V, B, Ry_hat_inv)
% Y: received data; V: carriers; B: watermarks
B_hat = sign(pinv(sqrtm(Ry_hat_inv)*V)*sqrtm(Ry_hat_inv)*Y);
ber_mat = (B_hat == B);
error_rate = length(find(ber_mat==0))/(size(B,1)*size(B,2));
end