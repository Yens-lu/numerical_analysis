import numpy as np
import pandas as pd
def SOR(A,b,omega,itera=1000):
    sor = [0,0]
    row = A.shape[0]
    x0 = np.zeros((row, 1))
    for i in range(itera):
        x1 = x0.copy()
        for k in range(row):
            a = A[k, 0:k].dot(x0[0:k]) + A[k, k + 1:row].dot(x0[k + 1:row])
            x0[k] = (1 - omega) * x0[k] + omega * (b[k, 0] - a) / A[k, k]
        theta = abs(x0[:] - x1[:])
        if theta.max() <= 0.000005:
            sor = [omega, i]
            break
        else :
            sor = [omega,i]
    return x0, sor
A = np.loadtxt('chapter3_matrix_R')
b = np.loadtxt('chapter3_matrix_V').reshape((A.shape[0],1))
S = pd.DataFrame(np.ones((2,2)),columns=['松弛因子','迭代次数'])
X = pd.DataFrame(np.ones((2,9)),columns=['I1 =','I2 =','I3 =','I4 =','I5 =','I6 =','I7 =','I8 =','I9 ='])
for n in range(1,100,1):
    omega = n/50
    re = SOR(A,b,omega)
    S.loc[n-1,:] = np.array(re[1]).reshape((1,2)).round(5)
    X.loc[n-1,:] = np.array(re[0]).T
print('解向量为：','\n',X.iloc[5,:].round(decimals=5),'\n')
S1 = S.sort_values(by='迭代次数')
print('最佳松弛因子','\n',S1.iloc[0,0],'\n'+'迭代次数为：'+'\n',S1.iloc[0,1])
#S.to_csv('iteration.txt',index='None',sep=',')