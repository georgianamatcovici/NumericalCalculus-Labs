#Cerinta 1: dimensiunea sistemului
import copy
import math
p_exp=int(input("p: "))
eps=10**-p_exp
solvable=True
p=0
q=0

def norma(x1, x2):
    res=0
    for i in range(len(x1)):
        res=res+(x1[i]-x2[i])**2
    return math.sqrt(res)
    
limit=10**10
k_max=10000

for i in range(1,5):
        d0=[]
        b=[]
        d1=[]
        d2=[]
        filename_d0=f"d0_{i}.txt"
        filename_b=f"b_{i}.txt"
        filename_d1=f"d1_{i}.txt"
        filename_d2=f"d2_{i}.txt"
        with open(filename_d0, 'r') as f:
            for line in f:
                if line.strip():
                  nr=float(line.strip())
                  if(nr<=eps):            #Cerinta 3
                      print("Zero pe diagonala principala sistemul ", i)
                      solvable=False
                      break
                  else : d0.append(nr)
            print(i, " d0 : ", len(d0))
        with open(filename_b, 'r') as f:
            for line in f:
                if line.strip():
                  b.append(float(line.strip()))
            n=len(b)
            print(i, " b : ", n)

        with open(filename_d1, 'r') as f:
            for line in f:
                if line.strip():
                  d1.append(float(line.strip()))
            p=n-len(d1)
            print(i, " p : ", p)

        with open(filename_d2, 'r') as f:
            for line in f:
                if line.strip():
                  d2.append(float(line.strip()))
            q=n-len(d2)
            print(i, " q : ", q)
        
        #cerinta 4
        xc=[0.0]*n
        xp=[0.0]*n
        k=0

        while True:
          for j in range(n):
            xp[j] = xc[j]

         
          for i1 in range(n):
            sum=0.0
            if (i1 - q >= 0): sum += d2[i1-q] * xc[i1-q]
            if (i1 - p >= 0): sum += d1[i1-p] * xc[i1-p]
            if (i1 + p < n):  sum += d1[i1] * xc[i1+p]
            if (i1 + q < n):  sum += d2[i1] * xc[i1+q]
            xc[i1] = (b[i1] - sum) / d0[i1]

          delta_x=norma(xc, xp)
          k=k+1
          if delta_x<eps:
             break
          if not (k<=k_max and delta_x<=limit):
           break
        if delta_x<eps:
            if(i==1):
             print("xc aproximare x*: ",xc)
            y = [0.0] * n
            norma=0.0
            for i in range(n):
                res_i = d0[i] * xc[i]
                if i - p>= 0:
                    res_i += d1[i-p]*xc[i-p]

                if i + p < n:
                    res_i += d1[i]*xc[i+p]
                 
                if i - q >= 0:
                    res_i += d2[i-q]*xc[i-q]
                
                if i + q < n:
                    res_i += d2[i] * xc[i+q]
                
               
                y[i] = res_i
                eroare=abs(y[i]-b[i])
                print("eroare: ", eroare)
        else: print("Divergenta")


