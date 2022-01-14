import numpy as np
import pandas as pd
#y(x)函数
def y(x):
    a = 3/(1+x**3)
    return a
#dy/dx = f(x,y)
def f(x,y):
    f = -(x**2)*y**2
    return f
#RK4通用程序
def RK4(x,y,h,i):
    k1 = f(x[i],y[i])
    k2 = f(x[i]+1/2*h,y[i]+1/2*h*k1)
    k3 = f(x[i]+1/2*h,y[i]+1/2*h*k2)
    k4 = f(x[i]+h,y[i]+h*k3)
    yi = y[i]+h/6*(k1+2*k2+2*k3+k4)
    y = np.append([y],[yi])
    return y
#AB4通用程序
def AB4(x,y,h,i):
    yi = y[i]+h/24*(55*f(x[i],y[i])-59*f(x[i-1],y[i-1])+37*f(x[i-2],y[i-2]) \
                      -9*f(x[i-3],y[i-3]))
    y = np.append([y],[yi])
    return y
#AB4-AM4预测校正通用程序
def AB4_AM4(x,y,h,i):
    yi0 = y[i]+h/24*(55*f(x[i],y[i])-59*f(x[i-1],y[i-1])+37*f(x[i-2],y[i-2]) \
                     -9*f(x[i-3],y[i-3]))
    yi = y[i]+h/24*(9*f(x[i+1],yi0)+19*f(x[i],y[i])-5*f(x[i-1],y[i-1])+f(x[i-2],y[i-2]))
    y = np.append([y],[yi])
    return y
#改进AB4-AM4预测校正通用程序
def speed_AB4AM4(x,y,h,i):
    yi0 = y[i] + h /24 * (55 * f(x[i], y[i]) - 59 * f(x[i - 1], y[i - 1]) + 37 * f(x[i - 2], y[i - 2])
                            - 9 * f(x[i - 3], y[i - 3]))
    yi1 = y[i] + h / 24 * (9 * f(x[i + 1], yi0) + 19 * f(x[i], y[i]) - 5 * f(x[i - 1], y[i - 1])
                           + f(x[i - 2], y[i - 2]))
    yi = 251/270*yi1 + 19/270*yi0
    y = np.append([y], [yi])
    return y

yx = np.round(np.array([3]),8)
x = np.round(np.array([j*0.1 for j in range(16)]),8)
result = pd.DataFrame(np.ones([x.shape[0],4]),index=['x'+str(i) for i in range(16)],
                      columns=['xi','yi','y(xi)','|y(xi)-yi|'])
result.iloc[0,1] = 3
result.iloc[:,0] = x
yrk = yx
#精确值
for i in range(16):
    result.iloc[i,2] = y(x=x[i])
#RK4的值
for i in range(0,15,1):
    m = i
    yrk = RK4(x,yrk,h=0.1,i=m)
result.iloc[:,1] = yrk
result.loc[:,'|y(xi)-yi|'] = abs(result.iloc[:,1]-result.iloc[:,2])
#AB4的值，AB4-AM4的值，speed-AB4AM4的值
yab = result.iloc[:4,1]
yabm = result.iloc[:4,1]
ysab = result.iloc[:4,1]
ab_result = result.copy()
abm_result = result.copy()
sab_result = result.copy()
for j in range(3,15):
    yab = AB4(x,yab,h=0.1,i=j)
    yabm = AB4_AM4(x,yabm,h=0.1,i=j)
    ysab = speed_AB4AM4(x,ysab,h=0.1,i=j)
ab_result.iloc[:,1] = yab
ab_result.loc[:,'|y(xi)-yi|'] = abs(ab_result.iloc[:,1]-ab_result.iloc[:,2])
abm_result.iloc[:,1] = yabm
abm_result.loc[:,'|y(xi)-yi|'] = abs(abm_result.iloc[:,1]-abm_result.iloc[:,2])
sab_result.iloc[:,1] = ysab
sab_result.loc[:,'|y(xi)-yi|'] = abs(sab_result.iloc[:,1]-sab_result.iloc[:,2])
result.to_csv('chapter6_result_of_RK4.csv',float_format='%.8f')
ab_result.to_csv('chapter6_result_of_AB4.csv',float_format='%.8f')
abm_result.to_csv('chapter6_result_of_AB4_AM4.csv',float_format='%.8f')
sab_result.to_csv('chapter6_result_of_speed_AB4AM4.csv',float_format='%.8f')


