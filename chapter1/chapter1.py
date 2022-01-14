'''第一章：舍入误差与有效数'''
import numpy as np
s = [pow(10,2),pow(10,4),pow(10,6)]
for n in s:
    N = n
    S = np.arange(2, N + 1, 1)
    S1 = np.square(S, dtype='float32') - 1
    S2 = 1 / S1
    SN = np.zeros([1, 2])
    SN.dtype = 'float32'
    '''由大到小'''
    for j in range(0, N - 1, 1):
        J = S2[j]
        SN[0, 0] = J + SN[0, 0]
    '''由小到大'''
    for j in range(0, N - 1, 1):
        J = S2[N - 2 - j]
        SN[0, 1] = J + SN[0, 1]
    print("当 N = %d 时" % N + '\n' + "由大到小的近似值: %.9f" % SN[0, 0] + '\n'
          + '有小到大的近似值: %.9f' % SN[0, 1],'\n')
#刚开始计算时没有将整数都设为单精度，导致最后的相加结果为负值