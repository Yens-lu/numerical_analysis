import numpy as np
import pandas as pd
def Gauss(A,b,p=1,s=1,I=range(10)):
    r = A.shape[0]
    Ab = np.hstack((A,b))
    x = np.zeros((r,1))
    for j in range(r):
        r0 = r - j - 1
        Ab[j:, j:] = Ab[j:, j:][np.argsort(-Ab[j:, j])]
        for i in range(r0):
            a = Ab[j + i + 1, j] / Ab[j, j]
            Ab[j + i + 1, j:] = Ab[j + i + 1, j:] - a * Ab[j, j:]
        A2 = Ab.copy()
    for n in range(r):
        k = r - 1 - n
        x[k, 0] = (Ab[k, r] - Ab[k, k + 1:r].sum(0)) / Ab[k, k]
        Ab[:k, k] = Ab[:k, k] * x[k, 0]
    Px = pd.DataFrame(x,index=I, columns=['各路电流为：']).round(decimals=5)
    if p == 1:
        print(Px, '\n')
    if s == 1:
        np.savetxt('Gauss.txt', A2, fmt='%.8f')
A = np.loadtxt('chapter3_matrix_R')
b = np.loadtxt('chapter3_matrix_V').reshape((A.shape[0],1))
Ix = ['I1 =','I2 =','I3 =','I4 =','I5 =','I6 =','I7 =','I8 =','I9 =']
Gauss(A,b,s=0,I=Ix)
