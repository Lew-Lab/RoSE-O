function out = proxmu2(x, recovStruct)

N = recovStruct.n_grid_p;
lambda = recovStruct.reg_val;
w = recovStruct.w_it;
n_f = recovStruct.subframe_l;

X = reshape(x, N, 12, n_f);
% X_s=squeeze(sqrt(sum(X.^2,2)))+eps;
X_s = reshape(sqrt(sum(X.^2, 2)), N, n_f) + eps;
w = repmat(w, N, 1);
q = max(X_s-w.*lambda, 0);
out = reshape(bsxfun(@times, reshape(q ./ X_s, N, 1, n_f), X), 12*N, n_f);
indx = out(1:3*N, :) <= 0;
out(indx > 0) = 0;
