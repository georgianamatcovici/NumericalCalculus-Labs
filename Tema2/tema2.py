import random

import numpy as np
import copy
from scipy.linalg import lu, lu_factor, lu_solve

flag=int(input("1-fisier 2-consola:")) #1-citesc din input.txt matricea din exemplul din tema  2-citesc n, eps de la tastatura, generez o matrice simetrica de n*n
if(flag==1):

    with open("input.txt", "r") as f:
        lines = f.readlines()

    n = int(lines[0].strip())
    eps = float(lines[1].strip())

    A = []
    for i in range(n):
        row = list(map(float, lines[2 + i].split()))
        A.append(row)
    b=[]
    b=list(map(float, lines[2+n].split()))

else:
    #generare matrice A simetrica si pozitiva de dimensiune n ca B*BT
    n = int(input("n:"))
    eps=float(input("eps:"))
    B = [[0]*n for i in range(0, n)]

    for i in range(0, n):
         for j in range(0, n):
             B[i][j]=random.random()

    A=[[0]*n for i in range(0, n)]

    for i in range(0, n):
         for j in range(0, n):
              s = 0
              for k in range(n):
                    s += B[i][k] * B[j][k]
              A[i][j] = s

    for i in range(n):
        A[i][i] += n
     

  #generare b random
    b = []

    for i in range(n):
        b.append(random.random())

    if(flag==1):
        print(b)

    A_copy = copy.deepcopy(A)

print("n =", n)
print("eps =", eps)

if(flag==1):
    print("A =\n", A)
    print("b=", b)


#Cerinta 1: Descompunere LU cu biblioteca scipy
P, L, U = lu(A)

if(flag==1):
    print("P =\n", P)
    print("L =\n", L)
    print("U =\n", U)


lu_fact, piv = lu_factor(A)
xlib = lu_solve((lu_fact, piv), b)
if(flag==1):
    print("xlib =", xlib)


#Cerinta2: descompunere LDLT cu Choleski
d=[0.0]*n

for p in range (0, n):
    sum=0.0
    #print(A[p][p])
    for k in range(0, p):
         sum=sum+d[k]*A[p][k]*A[p][k]
        #  print(A[k][p])
       

    # print(p)
    d[p]=A[p][p]-sum
    if abs(d[p]) < eps:
        print("Matricea nu este pozitiv definita.")
        break

    for i in range(p+1, n):
        sum=0.0
        for k in range(0, p):
            sum=sum+d[k]*A[i][k]*A[p][k]
        if(abs(d[p])>eps):
            A[i][p]=(A[i][p]-sum)/d[p]
        else:
            print("Nu poate fi calculata descompunerea")
            break
#Dupa aplicarea alg, in A in partea inferioara este retinuta partea inferioara a lui L, in d elementele de pe diagonala lui D



#Cerinta 3: determinantul atricei A care este produsul elementelor lui D

det=1

for i in range (0, n):
    det=det*d[i]
    
print("Determinantul: ", det)


 #Cerinta 4: calculul x_chol
z=[0.0]*n
for i in range(0, n):
    sum=0.0
    for j in range(0, i):
        sum=sum+A[i][j]*z[j]
    z[i]=b[i]-sum


y=[0.0]*n

for i in range(0, n):
    y[i]=z[i]/d[i]



x_chol=[0.0]*n
for i in range(n - 1, -1, -1):
    sum=0.0
    for j in range(i+1, n):
        sum=sum+A[j][i]*x_chol[j]
    x_chol[i]=y[i]-sum

if(flag==1):
    print("x_chol", x_chol)


#Cerinta 5: Calculul normelor

dif1=[0]*n
for i in range(0, n):
    res=0
    for j in range(0, n):
        if i<=j:
          res=res+A[i][j]*x_chol[j]
        else:
            res=res+A[j][i]*x_chol[j]
          
    dif1[i]=res-b[i]

    if(flag==1):
     print(dif1)

norma1=0

for i in range (0, n):
    norma1=norma1+dif1[i]**2
print("Norma1: ", np.sqrt(norma1))

norma2=0
for i in range (0, n):
    norma2=norma2+(xlib[i]-x_chol[i])**2
print("Norma2: ", np.sqrt(norma2))
