# Solving n x n sparse linear system of equations by LU Decomposition Iterative Refinement
# Assumptions: All diagonal entries are non-zeros
#              Tolerance for convergence is 0.001

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
    
# Getting the Constant vector
print ("Enter the constant array ->:",end=' ')
b=[int(x) for x in input().split()]

# Determining set P for fill-ins
P=[]
def get_fillins(B):
    for i in range(n):
        for j in range(n):
            if B[i][j]:
                P.append((i,j))

# Computing Boolean power up to level 1                
C=[[False]*n for i in range(n)]
for i in range(n):
    for j in range(n):
        for k in range(n):
            if B[i][k] and B[k][j]:
                C[i][j]=True
                break
get_fillins(C)

# Updating LUa,LUi,LUj
LUa,LUi,LUj=VA[:],IA[:],JA[:]           # The product (L*U)
if P!=P_temp:
    d={}
    LUi=[0]*(n+1)
    for i in range(n):d[i]=0
    for k,(i,j) in enumerate(P):
        d[i]+=1
        if len(P_temp)<=k or P_temp[k]!=P[k]:
            LUa.insert(k,0)
            LUj.insert(k,j)
            P_temp.insert(k,(i,j))
    for i in range(n):
        LUi[i+1]=LUi[i]+d[i]

## Diagonal entry locations of U
Diag=[]
for i in range(n):
    for j in range(LUi[i],LUi[i+1]):
        if LUj[j]==i:Diag.append(j)

## Incomplete LU Decomposition       
Point=[0]*n
for i in range(1,n):
    for v in range(LUi[i]+1,LUi[i+1]):Point[LUj[v]]=v
    for v in range(LUi[i],Diag[i]):
        j=LUj[v]
        LUa[v]=LUa[v]/LUa[Diag[j]]
        for w in range(Diag[j]+1,LUi[j+1]):
            k=Point[LUj[w]]
            if k>0:
                LUa[k]=LUa[k]-LUa[v]*LUa[w]            
    for v in range(LUi[i]+1,LUi[i+1]):
        Point[LUj[v]]=0

def Solve(x,b):            	              # For any (LU)*x=b, Solving for x
    ## Forward substitution Ly=b
    y=[0 for i in range(n)]
    y[0]=b[0]
    for i in range(1,n):
        for j in range(LUi[i],Diag[i]):
            y[i] += LUa[j]*y[LUj[j]]
        y[i]=b[i]-y[i]
    
    ## Back substitution Ux=y
    x[n-1]=y[n-1]/LUa[Diag[n-1]]
    for i in range(n-2,-1,-1):
        for j in range(Diag[i],LUi[i+1]):
            x[i] += LUa[j]*x[LUj[j]]
        x[i]=(y[i]-x[i])/LUa[Diag[i]]
x=[0 for i in range(n)] 
Solve(x,b)     

def is_converged(x_temp,x):
    tolerance=0.001
    for i in range(n):
        if abs(x_temp[i]-x[i])>tolerance:return False
    return True

x_temp=[0]*n
iterations=0
while (not is_converged(x_temp,x)) and iterations<=1000:
    iterations+=1
    x_temp=x[:]

    ## Compute r = b - A*x
    Ax=[0]*n
    for i in range(n):
        tot=0
        for v in range(IA[i],IA[i+1]):
            tot+=VA[v]*x_temp[JA[v]]
        Ax[i]=tot        
    r=[0]*n
    for i in range(n):r[i]=b[i]-Ax[i]

    ## (LU)*d = r, Solving for d
    d=[0]*n
    Solve(d,r)

    ## Compute x2 = x1 + d
    for i in range(n):x[i]+=d[i]

if iterations<=1000:print("The solution vector X is",x)
else:print("The solution diverges by this method")
            
