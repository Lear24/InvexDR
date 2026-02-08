function  S = util_CollVecBeta(beta_cell,m,p,d)

S = zeros(p*d,m);

for i =1:m
    temp = cell2mat(beta_cell(i));
    S(:,i) = temp(:);
end


end