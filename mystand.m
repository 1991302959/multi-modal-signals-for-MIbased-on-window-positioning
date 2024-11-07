 function out=mystand(A)

out = [];
n=size(A);                  %获取矩阵大小
minA = min(A,[],'all');     %获取矩阵最小值
maxA = max(A,[],'all');     %获取矩阵最大值
out = (A - repmat(minA,n))./(maxA-minA);    %最大最小值归一化

