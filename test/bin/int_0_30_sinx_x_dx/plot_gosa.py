import matplotlib.pyplot as plt
import numpy as np

s=[]
g10=[]
f10="./fort.10"
data = np.loadtxt(f10, comments="#")
s=data[:,0]
g10=data[:,1]

g11=[]
f11="./fort.11"
data = np.loadtxt(f11, comments="#")
g11=data[:,1]

g12=[]
f12="./fort.12"
data = np.loadtxt(f12, comments="#")
g12=data[:,1]

g13=[]
f13="./fort.13"
data = np.loadtxt(f13, comments="#")
g13=data[:,1]

g14=[]
f14="./fort.14"
data = np.loadtxt(f14, comments="#")
g14=data[:,1]

g15=[]
f15="./fort.15"
data = np.loadtxt(f15, comments="#")
g15=data[:,1]

g16=[]
f16="./fort.16"
data = np.loadtxt(f16, comments="#")
g16=data[:,1]

g17=[]
f17="./fort.17"
data = np.loadtxt(f17, comments="#")
g17=data[:,1]

g18=[]
f18="./fort.18"
data = np.loadtxt(f18, comments="#")
g18=data[:,1]

g19=[]
f19="./fort.19"
data = np.loadtxt(f19, comments="#")
g19=data[:,1]

# g20=[]
# f20="./fort.20"
# data = np.loadtxt(f20, comments="#")
# g20=data[:,1]

plt.xlabel('The number of nodes (n)') # x軸のラベル
plt.ylabel('Relative $\ell_2$ error') # y軸のラベル

plt.title('$\int_{0}^{30} \sin(x)/x \mathrm{d}x$')
plt.yscale('log')
plt.ylim(10**(-16), 10**6)
plt.grid(linestyle='--')
plt.plot(s, g10, marker='o', linestyle='-', label='$s=0$')
plt.plot(s, g11, marker='o', linestyle='-', label='$s=1$')
plt.plot(s, g12, marker='o', linestyle='-', label='$s=2$')
plt.plot(s, g13, marker='o', linestyle='-', label='$s=3$')
plt.plot(s, g14, marker='o', linestyle='-', label='$s=4$')
plt.plot(s, g15, marker='o', linestyle='-', label='$s=5$')
plt.plot(s, g16, marker='o', linestyle='-', label='$s=6$')
plt.plot(s, g17, marker='o', linestyle='-', label='$s=7$')
plt.plot(s, g18, marker='o', linestyle='-', label='$s=8$')
plt.plot(s, g19, marker='o', linestyle='-', label='$s=9$')
plt.legend()
plt.savefig('int_0_30_sinx_dx_x.eps')
