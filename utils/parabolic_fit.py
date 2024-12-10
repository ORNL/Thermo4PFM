import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#usage example:
#python parabolic_fit.py phaseL 900. FvsC900.csv 920. FvsC920.csv 0.7 0.85

#
def poly2(x, A, B, C):
  y = A*x*x+B*x+C
  return y

#############################

nargs = len(sys.argv)

phase = sys.argv[1]
print("phase {}".format(phase))

if phase=='phaseL':
  col = 1
if phase=='phaseA':
  col = 2
if phase=='phaseB':
  col = 3

filenames = []
temperatures = []

temperatures.append(sys.argv[2])
temperatures.append(sys.argv[4])

Tref = eval(temperatures[0])

filenames.append(sys.argv[3])
filenames.append(sys.argv[5])

a=[]
b=[]
c=[]

for filename, temperature in zip(filenames,temperatures):
  print("filename = {}".format(filename))
  print("temperature = {}".format(temperature))
  ra = np.genfromtxt(filename, delimiter=',',dtype=None, names=True)

  aa = ra.view(np.float64).reshape(len(ra), -1)

  nrows = aa.shape[0]

  #limit x range on right
  if nargs>7:
    xmax = eval(sys.argv[7])
    print("max. value for x = {}".format(xmax))
    n=int(nrows*xmax)
    print(n)
    for i in range(nrows-n):
      #print("remove row {}".format(nrows-i-1))
      aa = np.delete(aa, (nrows-i-1), axis=0)

  #limit x range on left
  if nargs>6:
    xmin = eval(sys.argv[6])
    print("Min. value for x = {}".format(xmin))
    n=int(nrows*xmin)
    for i in range(n):
      #print('remove row 0')
      aa = np.delete(aa, (0), axis=0)

  x=aa[:,0]

  y=aa[:,col]

  p, covariance = curve_fit(poly2, x,y)
  a.append(p[0])
  b.append(p[1])
  c.append(p[2])
  print("T = {}, polynomial {}*x^2 + {}*x + {}".format(temperature,p[0],p[1],p[2]))

  fit_y = poly2(x, p[0], p[1], p[2])
  plt.plot(x,y,'o')
  plt.plot(x,fit_y,'-')
  #plotfile = "parabolic_fit"+str(temperature)+".png"
  #plt.savefig(plotfile, dpi=100)

#get temperature dependent coefficients
dTinv=1./(eval(temperatures[1])-eval(temperatures[0]))
a0=2*a[0]
a1=2*(a[1]-a[0])*dTinv
b0=b[0]
b1=(b[1]-b[0])*dTinv
c0=c[0]
c1=(c[1]-c[0])*dTinv
print("0.5*({}+{}*(T-{}))*x^2 + ({}+{}*(T-{}))*x + ({}+{}*(T-{}))".format(a0,a1,Tref,b0,b1,Tref,c0,c1,Tref))

#check result for intermediate temperature
temperature=0.5*(eval(temperatures[0])+eval(temperatures[1]))
p0=0.5*(a0+a1*(temperature-Tref))
p1=b0+b1*(temperature-Tref)
p2=c0+c1*(temperature-Tref)

fit_y = poly2(x, p0, p1, p2)
plt.plot(x,fit_y,'-')
plotfile = "parabolic_fit"+str(temperature)+".png"
plt.savefig(plotfile, dpi=100)

