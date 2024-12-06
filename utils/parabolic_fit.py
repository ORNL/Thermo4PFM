import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#
def poly2(x, A, B, C):
  y = A*x*x+B*x+C
  return y

nargs = len(sys.argv)

ra = np.genfromtxt('FvsC.csv', delimiter=',',dtype=None, names=True)

aa = ra.view(np.float64).reshape(len(ra), -1)

nrows = aa.shape[0]

#limit x range on left
if nargs>1:
  xmin = eval(sys.argv[1])
  print("Min. value for x = {}".format(xmin))
  n=int(nrows*xmin)
  for i in range(n):
    print('remove row 0')
    aa = np.delete(aa, (0), axis=0)

#limit x range on right
if nargs>2:
  xmax = eval(sys.argv[2])
  print("max. value for x = {}".format(xmax))
  n=int(nrows*xmax)
  print(n)
  for i in range(nrows-n):
    print("remove row {}".format(nrows-i-1))
    aa = np.delete(aa, (nrows-i-1), axis=0)


x=aa[:,0]

#loop over columns > 0 (phases)
for j in range(1,aa.shape[1]):
  y=aa[:,j]

  parameters, covariance = curve_fit(poly2, x,y)
  fit_A=parameters[0]
  fit_B=parameters[1]
  fit_C=parameters[2]
  print("polynomial {}*x^2 + {}*x + {}".format(fit_A,fit_B,fit_C))

  fit_y = poly2(x, fit_A, fit_B, fit_C)
  plt.plot(x,y,'o')
  plt.plot(x,fit_y,'-')

plt.savefig('parabolic_fit.png', dpi=100)

