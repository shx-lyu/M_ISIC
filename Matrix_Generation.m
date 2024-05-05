function Mat_Gen = Matrix_Generation
Mat_Gen.Gen_U = @Gen_U;
Mat_Gen.Gen_B = @Gen_B;
end

% N: Number of generated U
% sourceNum and len denote the width and height of U, respectively
function Uc = Gen_U(N,sourceNum,len)
Uc = cell(1,N);
for i = 1:N
    U = zeros(len,sourceNum);
    for j = 1:sourceNum
        U(:,j) = randn(len,1);
        U(:,j) = U(:,j)/norm(U(:,j));
    end
    Uc{1,i} = U;
end
end

function Bc = Gen_B(P,sourceNum,w_len)
Bc = cell(1,P);
for i = 1:P
    B_tmp = zeros(sourceNum,w_len);
    for j = 1:sourceNum
        B_tmp(j,:) = sign(randn(1,w_len));
    end
    Bc{1,i} = B_tmp;
end
end
