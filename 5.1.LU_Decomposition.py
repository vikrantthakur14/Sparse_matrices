# Solving n x n sparse linear system of equations by LU Decomposition
# Assumptions: The coefficient matrix is arranged such that all diagonal elements are non-zero

print ("The number of variables in the system :",end=' ')
n=int(input())
print("Enter the (%d) x (%d) coefficient matrix :"%(n,n))
NNZ=0
VA,JA=[],[]
IA=[0]*(n+1)
B=[[False]*n for i in range(n)]         # Boolean matrix
P_temp=[]
for i in range(n):
    temp=0
    for j,value in enumerate(map(float,input().split())):
        if value!=0:
            temp+=1
            VA.append(float(value))
            JA.append(j)
            B[i][j]=True
            P_temp.append((i,j))
    IA[i+1] = IA[i] + temp
    NNZ+=temp
    
## Getting the Constant vector
print ("Enter the constant array ->:",end=' ')
b=[int(x) for x in input().split()]

## Determining set P for fill-ins
P=[]
def get_fillins(B):
    for i in range(n):
        for j in range(n):
            if B[i][j]:
                P.append((i,j))

## Applying Warshall's Algorithm 
for k in range(n):
    for i in range(n):
        for j in range(n):
            B[i][j]=(B[i][j] or (B[i][k] and B[k][j]))
get_fillins(B)

## Updating VA,IA,JA
if P!=P_temp:
    d={}
    IA=[0]*(n+1)
    for i in range(n):d[i]=0
    for k,(i,j) in enumerate(P):
        d[i]+=1
        if len(P_temp)<=k or P_temp[k]!=P[k]:
            VA.insert(k,0)
            JA.insert(k,j)
            P_temp.insert(k,(i,j))
    for i in range(n):
        IA[i+1]=IA[i]+d[i]

## Diagonal entry locations of U
Diag=[]
for i in range(n):
    for j in range(IA[i],IA[i+1]):
        if JA[j]==i:Diag.append(j)

## Incomplete LU Decomposition        
Point=[0]*n
for i in range(1,n):
    for v in range(IA[i]+1,IA[i+1]):Point[JA[v]]=v
    for v in range(IA[i],Diag[i]):
        j=JA[v]
        VA[v]=VA[v]/VA[Diag[j]]
        for w in range(Diag[j]+1,IA[j+1]):
            k=Point[JA[w]]
            if k>0:
                VA[k]=VA[k]-VA[v]*VA[w]            
    for v in range(IA[i]+1,IA[i+1]):
        Point[JA[v]]=0
        
## Forward substitution Ly=b
y=[0 for i in range(n)]
y[0]=b[0]
for i in range(1,n):
    for j in range(IA[i],Diag[i]):
        y[i]=y[i]+VA[j]*y[JA[j]]
    y[i]=b[i]-y[i]

## Back substitution Ux=y
x=[0 for i in range(n)]
x[n-1]=y[n-1]/VA[Diag[n-1]]
for i in range(n-2,-1,-1):
    for j in range(Diag[i],IA[i+1]):
        x[i]=x[i]+VA[j]*x[JA[j]]
    x[i]=(y[i]-x[i])/VA[Diag[i]]    

print("The solution vector X is",x)
