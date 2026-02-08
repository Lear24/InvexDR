function X_cell = util_vector(X_cell,m,n,p,d)

X_cell_temp = cell(m,1);

for i = 1:m
    X = cell2mat(X_cell(i,1));
    X1 = zeros(p*d,n(i));
    temp = zeros(p,d);
    for j = 1:n(i)
       temp = X(:,(j-1)*d+1:j*d);
       X1(:,j)=temp(:);
    end
    X_cell_temp(i,1) = mat2cell(X1,size(X1,1),size(X1,2));
    
end
X_cell = X_cell_temp;

end