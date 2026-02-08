function X_cell = util_invvector(X_cell,m,n,p,d)

X_cell_temp = cell(m,1);

for i = 1:m
    X = cell2mat(X_cell(i,1)); % p*d-by-n
    X1 = zeros(p,n(i)*d);
    for j = 1:n(i)
       X1(:,(j-1)*d+1:j*d) = reshape(X(:,j),p,d);
    end
    X_cell_temp(i,1) = mat2cell(X1,p,n(i)*d);
end
X_cell = X_cell_temp;

end