function c = FDRcutoff(p,alpha,flag)
% function cutoff = FDRcutoff(pValues,alpha,flag=false)
% 
% pValues is a vector of valid p-values
% 
% rejecting those hypotheses corresponding to pValues <= cutoff(k)
% controls the FDR at level alpha(k)
% 
% setting flag=true assumes that the pValues are independent (or more
% generally, positively dependent in a certain sense as specified in the
% FDR literature)

if nargin < 3 || isempty(flag), flag = false; end

p = sort(p(:));
D = numel(p);

if flag, alpha = alpha/D; else alpha = alpha/(D*sum(1./(1:D))); end

ndx = bsxfun(@le,p./(1:D).',alpha(:).');
a = any(ndx,1);

c = zeros(size(alpha));

for k = 1:numel(alpha)
    if a(k)
        c(k) = p(find(ndx(:,k),1,'last'));
    end
end