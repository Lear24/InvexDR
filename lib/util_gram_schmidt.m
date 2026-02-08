function Q=util_gram_schmidt(A)
% 实现 Schmidt 正交化
% 向量以列向量的形式进行排列, 因此矩阵维度为 d,n

[p,d]=size(A);
Q = zeros(p,d);
R = zeros(p,p);

for j=1:d
    v = A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        v = v - R(i,j)'*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j) = v/R(j,j);
end


end