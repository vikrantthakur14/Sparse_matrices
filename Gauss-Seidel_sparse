# Finding Solution of a n x n sparse linear system of equations using Gauss-Seidel method
# Compressed Row Storage is used for storing the sparse coefficient matrix
# Assumption -> Given coefficient matrix is diagonally dominant so as the method converges
#               Tolerance = 0.001 and if iterations exceed 1000 No feasible solution (i.e method diverges)

# CRS for coeff matrix
print("The number of variables in the system :",end=' ')
n=int(input())
print("Enter the %d x %d coefficient matrix which is diagonally dominant:"%(n,n))
NNZ=0
VA,JA=[],[]
IA=[0]*(n+1)
for i in range(n):
    temp=0
    for j,value in enumerate(map(float,input().split())):
        if value!=0:
            temp+=1
            VA.append(value)
            JA.append(j)
    IA[i+1] = IA[i] + temp
    NNZ+=temp
    
# Getting the Constant matrix
print("Enter the constant array -> :",end=' ')
B=[float(x) for x in input().split()]

# Stopping the iterations
tolerance=1e-04
def is_converging(arr1,arr2):                             # To check if prev soln is same as new_soln
    for a,b in zip(arr1,arr2):
        if abs(abs(a)-abs(b))>tolerance:return False
    return True

# Gauss-Seidel algorithm
soln=[0]*n
new_soln=[1]*n
iters=0            
while not is_converging(soln,new_soln) or iters > 1000:   # Breaks when solution converges or number of iterations exceed 10^5
    iters+=1
    new_soln=soln[:]
    for x in range(n):
        tot=0
        for v in range(IA[x],IA[x+1]):
            if JA[v]!=x: tot += VA[v]*soln[JA[v]]
            else: at_diag=VA[v]            
        soln[x]=(B[x]-tot)/at_diag

if iters<=1000: print("The solution to the system of linear equations is :",soln)
else: print("The given system of equation diverges and therefore no feasible solution by this method")

        
    
    
    
