# CSR (Compressed sparse row) storage scheme for sparse square matrices
# Storing a given sparse (m x m) square matrix in CSR format
# Example: A = [[0 4 0 6]  ;(4 x 4) sparse square matrix
#               [1 8 0 0]   VA = [4,6,1,8,7,8,9]
#    NNZ=7      [0 0 7 0]   IA = [0,2,4,5,7]
#               [8 0 9 0]]  JA = [1,3,0,1,2,0,2]

m=int(input("Enter the order 'm' of your square matrix: "))

NNZ=0        # Number of non zero elements in matrix 'A'    
VA=[]        # array storing all non zero elements of matrix 'A'(length=NNZ)
IA=[0]*(m+1) # First m elements are the 'VA' indices of first non-zero occurences in each row
             # and the m+1 element represents the 'NNZ' (length=m+1)
JA=[]        # Column index in 'A' of each element of 'VA' (length=NNZ)

print("Input the %d x %d matrix: "%(m,m))
for i in range(m):
    temp=0
    for j,value in enumerate(map(float,input().split())):
        if value!=0:
            temp+=1
            VA.append(value)
            JA.append(j)
    IA[i+1] = IA[i] + temp
    NNZ+=temp

print("NNZ =",NNZ)
print("VA =",VA)
print("IA =",IA)
print("JA =",JA)
