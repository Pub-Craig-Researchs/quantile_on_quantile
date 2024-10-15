function loss = qqrLossFunc(beta,dataTbl,theta,tao,h,IF_LAGGED_IDV)
% dataTbl should be a table, in which the denpendent varibale should be
% allocated in the last column
dpV = dataTbl{:,end};
idVs = dataTbl{:,1:end-1};
alpha = beta(end);

if(IF_LAGGED_IDV)
    idVs = [NaN;idVs(1:end-1)];
end
lag_dpV = [NaN;dpV(1:end-1)];
regTbl = rmmissing(table(idVs,lag_dpV,dpV,VariableNames=...
    {'independent_Varibales','lagged_dependent_Variable','dependent_Variable'}));

fit = beta(1) + beta(2) * (regTbl.independent_Varibales - quantile(regTbl.independent_Varibales,tao)) + ...
    alpha * regTbl.lagged_dependent_Variable;

loss = regTbl.dependent_Variable - fit;
loss(loss>0) = loss(loss>0) * theta;
loss(loss<0) = loss(loss<0) * (theta - 1);

Idx = tiedrank(regTbl.independent_Varibales,1);
F = (Idx-1) / height(regTbl);
K = normpdf((F - tao)/h);
loss = sum(loss .* K);
end
