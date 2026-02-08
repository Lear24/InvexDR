function beta=set_initialization(m,p,d,varargin)
% varargin = 'binomial' or  'normal'
% beta is a m-by-1 cell with p-by-d matrix
mode = varargin{1};
beta = cell(m,1);

switch mode
    case 'binomial'
        for i = 1:m
            alpha = binornd(1,0.5,p,d);
            beta(i)=mat2cell(alpha,p,d);
        end
    case 'normal'
        for i = 1:m
            alpha = randn(p,d);
            beta(i)=mat2cell(alpha,p,d);
        end
end

if ~isempty(varargin{2}); modeB = varargin{2};else; modeB = 'nope';end
switch modeB
    case 'orth'
        for i = 1:m
            alpha = beta{i};
            [~,J]=rref(alpha); % J 中为A的极大线性无关组的坐标信息
            b = util_gram_schmidt(alpha(:,J)); % 正交化
            beta{i} = b;
        end
end



end