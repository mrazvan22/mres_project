function normMat = normaliseVarMatrix(mat)

nrBiomk = size(mat,1);
sumMat = sum(mat,2)
normMat = mat ./ repmat(sumMat, [1, nrBiomk]);

end