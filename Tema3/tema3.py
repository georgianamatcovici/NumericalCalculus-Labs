import math
import copy
import random
import numpy as np

x=int(input("1-fisier 2-consola:"))

if(x==1):


    with open("input.txt", "r") as f:
            lines = f.readlines()

    n = int(lines[0].strip())
    eps = float(lines[1].strip())

    A = []
    for i in range(n):
        row = list(map(float, lines[2 + i].split()))
        A.append(row)
        s=[]
        s=list(map(float, lines[2+n].split()))
else:
     n = int(input("n:"))
     eps=float(input("eps:"))
     
     A = [[0]*n for i in range(0, n)]
     for i in range(0, n):
         for j in range(0, n):
             A[i][j]=random.uniform(-100, 100)
     
     s=[0]*n

     for i in range(n):
         s[i]=random.random()

         
A_copy=np.array(A)
print("n =", n)
print("eps =", eps)
# print("A =\n", A)
# print("s=", s)

#Cerinta 1 - calcul b
b = []
for i in range(n):
    suma_rand = 0
    for j in range(n):
        suma_rand += A[i][j] * s[j]
    b.append(suma_rand)

if(x==1): print("b =", b)


b_copy=copy.deepcopy(b)


#Cerinta 2: descompunere QR folosind Householder
Q_T=[[0]*n for i in range(0, n)]
for i in range(0, n):
    Q_T[i][i]=1

u=[]
for r in range(0, n-1):
    sigma=0.0
   
    for i in range(r, n):
        sigma=sigma+abs(A[i][r])**2
    # print(sigma)
    if(sigma<=eps):
        break
    k=math.sqrt(sigma)
    if(A[r][r].real>0):
        k=-k
    beta=sigma-k*A[r][r]
    u=[0.0]*n
    u[r]=A[r][r]-k
    for i in range(r+1, n):
        u[i]=A[i][r]

    for j in range(r+1, n):
        gama=0
        for i in range(r, n):
            gama=gama+u[i]*A[i][j]
        gama/=beta

        for i in range(r, n):
            A[i][j]=A[i][j]-gama*u[i]

    A[r][r]=k
    for i in range(r+1, n):
        A[i][r]=0
    
    gama=0

    for i in range(r, n):
        gama=gama+u[i]*b[i]
    gama=gama/beta

    for i in range(r, n):
        b[i]=b[i]-gama*u[i]
    
    for j in range(0, n):
        gama=0
        for i in range(r, n):
            gama=gama+u[i]*Q_T[i][j]
        gama=gama/beta

        for i in range(r, n):
            Q_T[i][j]=Q_T[i][j]-gama*u[i]
if(x==1):
    print("R",A)
    print("Q transpus", Q_T)


# R- superior triunghiulara stocata in A, Q_T transpusa lui Q, b=Q_T*b_init


# Cerinta 3: rezolvare sistem 

Q1, R1 = np.linalg.qr(A_copy)

if(x==1):
    print("Q1", Q1)
    print("R1", R1)


b_sistem=[0]*n
b_sistemQR=[0]*n


for i in range(0, n):
    sum=0 
    sum1=0
    for j in range(0, n):
        sum=sum+Q_T[i][j]*b_copy[j]
        sum1=sum1+Q1[j][i]*b_copy[j]
    b_sistem[i]=sum
    b_sistemQR[i]=sum1


if(x==1):
    print("b", b_sistem)
    print("bQR", b_sistemQR)


x_householder=[0]*n
x_QR=[0]*n
for i in range(n-1, -1, -1):
    sum=0
    sum1=0
    for j in range(i+1, n):
        sum=sum+A[i][j]*x_householder[j]
        sum1=sum1+R1[i][j]*x_QR[j]
    x_householder[i]=(b_sistem[i]-sum)/A[i][i]
    x_QR[i]=(b_sistemQR[i]-sum1)/R1[i][i]

if(x==1):
    print("x", x_householder)
    print("x", x_QR)

norma=0
for i in range(0, n):
    norma=norma+(x_QR[i]-x_householder[i])**2

print("Norma solutiei x", norma)


#CErinta 4: erori
err1=0
err2=0
dif1=[0]*n
dif2=[0]*n

for i in range(0, n):
    res1=0
    res2=0
    for j in range(0, n):
          res1=res1+A_copy[i][j]*x_householder[j] 
          res2=res2+A_copy[i][j]*x_QR[j] 

    dif1[i]=res1-b_copy[i]
    dif2[i]=res2-b_copy[i]

    err1=err1+dif1[i]**2
    err2=err2+dif2[i]**2

print("Eroare1: ", math.sqrt(err1))
print("Eroare2: ", math.sqrt(err2))


err3=0

norma1=0
norma2=0
norma3=0
err4=0

for i in range(0, n):
    norma1+=(x_householder[i]-s[i])**2
    norma2+=(s[i])**2
    norma3+=(x_QR[i]-s[i])**2
  

err3=math.sqrt(norma1)/math.sqrt(norma2)
err4=math.sqrt(norma3)/math.sqrt(norma2)

print("Eroare3: ", err3)
print("Eroare4: ", err4)


#Cerinta 5: inversa

det=1
for i in range(0, n):
    if(abs(A[i][i])<eps):
         print("Nu putem calcula inversa")
         break
    det=det*A[i][i]

print("Determinant:", det)


Inv_A=[[0]*n for  i in range(0, n)]


for j in range(n):
    coloana_j = []
    x_star=[0.0]*n
    for i in range(n):
        coloana_j.append(Q_T[i][j])
    #print(coloana_j)
    for k in range(n-1, -1, -1):
        sum=0
        for l in range(k+1, n):
            sum=sum+A[k][l]*x_star[l]
        x_star[k]=(coloana_j[k]-sum)/A[k][k]
        Inv_A[k][j]=x_star[k]

if(x==1): 
 print("Inversa", Inv_A)

Inv_A_bibl = np.linalg.inv(A_copy)
if(x==1):
  print("Inversa Librarie", Inv_A_bibl)

norma_inversa=0
for i  in range(0, n):
    for j in range(0, n):
        norma_inversa+=(Inv_A[i][j]-Inv_A_bibl[i][j])**2
print("Norma intre inverse matrici", math.sqrt(norma_inversa))

