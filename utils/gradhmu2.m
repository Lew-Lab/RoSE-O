function out = gradhmu2(x, recovStruct)
mu = recovStruct.mu;
out = (x - proxmu2(x, recovStruct)) * (1 / mu);
end
