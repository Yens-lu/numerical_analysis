import numpy as np
import pandas as pd
def CSI(data,dy0=0.2,dyn=0.8): #求解插值函数的各个系数
    n = data.shape[1]
    result = pd.DataFrame(np.zeros([n, 3]), columns=['Mj', 'f[xj,xj+1]', 'hj'])
    D = np.round(np.zeros([1, n]), 6)
    b = np.round(np.ones([1, n]) * 2, 6)
    mu = np.round(np.ones([1, n - 1]), 6)
    lam = np.round(np.ones([1, n - 1]), 6)
    M = np.round(np.zeros([1, n]), 6)
    '''差商计算,用矩阵加减来算差商'''
    deltah = data[1, 1:] - data[1, :-1]
    f = data[2, 1:] - data[2, :-1]
    ff = f.reshape([1, n - 1]) / deltah.reshape([1, n - 1])
    D[0, 1:n - 1] = (ff[0, 1:] - ff[0, :-1]).reshape([1, n - 2]) / (data[1, 2:] - data[1, :-2]).reshape([1, n - 2])
    D[0, 0] = ((data[2, 1] - data[2, 0]) / (data[1, 1] - data[1, 0]) - dy0) / (data[1, 1] - data[1, 0])
    D[0, n - 1] = (dyn - (data[2, n - 1] - data[2, -2]) / (data[1, n - 1] - data[1, -2])) / (
                data[1, n - 1] - data[1, -2])
    mu[0, :-1] = (deltah[:-1] / (deltah[:-1] + deltah[1:])).reshape([1, n - 2])
    lam[0, 1:] = 1 - mu[0, :-1]
    '''追赶法（GAUSS消去）'''
    beta = b
    y = D
    for j in range(1, n):
        l = mu[0, j - 1] / beta[0, j - 1]
        beta[0, j] = beta[0, j - 1] - l * lam[0, j - 1]
        y[0, j] = y[0, j] - l * y[0, j - 1]
    M[0, n - 1] = y[0, n - 1] / beta[0, n - 1]
    for j in range(1, n):
        M[0, n - 1 - j] = (y[0, n - 1 - j] - lam[0, n - 1 - j] * M[0, n - j]) / beta[0, n - 1 - j]
    result.iloc[:,0] = M.T
    result.iloc[:-1, 2] = deltah[:].T
    result.iloc[:-1, 1] =  ff[0,:].T
    return result
def sx(x,co,data): #求解插值函数S(x)
    for i in range(10):
        if (x > data[1, i]) & (x <= data[1, i + 1]):
            S = data[2, i] + (co[i,1] - (1 / 3 * co[i,0] + 1 / 6 * co[i + 1,0]) * co[i,2]) * (x - data[1, i]) + \
                1 / 2 * co[i,0] * ((x - data[1, i]) ** 2) + 1 / (6 * co[i,2]) * (co[i + 1,0] - co[i,0]) * \
                ((x - data[1, i]) ** 3)
    return S
data = np.loadtxt('chapter4_data',delimiter=',')
#求解Mj,hj,f[xj,xj+1]
result = CSI(data)
coefficiency = np.array(result.iloc[:,:]).reshape([11,3])
indexs = []
for x in range(0,10,):
    indexs.append('S('+str(x+0.5)+')=')
#求解插值函数S(x)
Sxi = pd.DataFrame(np.zeros([10,1]),index=indexs,columns=['插值'])
for i in range(0,10):
    Sxi.iloc[i,0] = sx(x=i+0.5,co=coefficiency,data=data)
Sxi.to_csv('Chapter4_result',sep=' ',header=False,float_format='%.8f')
