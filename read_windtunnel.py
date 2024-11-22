import numpy as np
import matplotlib.pyplot as plt

file_path= r'C:\Users\frede\Downloads\raw_testG9.txt'
Line= 10

List=np.genfromtxt(file_path,delimiter='')
pos=[]
p_top=[]
p_bottom=[]
pos2=[]
x=32
b=0
for i in range(7,x):
    p_top.append(List[Line,i])
    pos.append(b)
    b+=0.667
b=0
for i in range(x,57):
    pos2.append(b)
    p_bottom.append(List[Line,i])
    b+=0.667

print(List[Line,2])
plt.plot(pos,p_top)
plt.plot(pos2,p_bottom)
plt.gca().invert_yaxis()
plt.grid()

plt.show()
