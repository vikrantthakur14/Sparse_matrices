# Finding Solution of a n variable sparse linear system of equations using Cramers rule
# Assumptions: Coefficient matrix is non-singular (i.e Determinant exists)

print("The number of variables in the system Ax=b:",end=' ')
n=int(input())
print("Enter the %d x %d Augmented matrix [A|b]:"%(n,n+1))
NNZ=0
VA={};JA=[]
IA=[0]*(n+1)

for i in range(n):
    temp=0
    for j,value in enumerate(map(float,input().split())):
        if value!=0:
            temp+=1
            VA[(i,j)]=value
            JA.append(j)
    IA[i+1] = IA[i] + temp
    NNZ+=temp

S=[]                   # Numerical structure of Augmented matrix [A|b]
for x in range(n):
    S.append(JA[IA[x]:IA[x+1]])

sorted_S=list()        # The rows of S are sorted in increasing length of set H
sorted_rows=list()     # row indices after sorting

for row_ind,H in sorted(enumerate(S),key=lambda x:len(x[1])):   
    sorted_S.append(H)                         # H -> set of column indices of non zero entries in each row in S
    sorted_rows.append(row_ind)    

# Taking permutations of each H such that column indices are unique    
unique_cols=[]
for i in sorted_S[0]:unique_cols.append([i])

for h in range(1,n):
    temp=len(unique_cols)
    for y in range(len(unique_cols)):
        for x in sorted_S[h]:
            if x not in unique_cols[y]:
                unique_cols.append(unique_cols[y]+[x])
    unique_cols=list(unique_cols[temp:])

# Creating numerical substructures
subs_S=[[] for x in range(n+1)]                # All substructures 0,1,...n
for x in range(n+1):
    for cols in unique_cols:
        if x not in cols:
            subs_S[x].append(cols)

# Finding determinants of all substructures 0,1...n
Det_subs=[]
def no_of_ways(arr1,arr3):                     # No. of interchanges needed in h
    arr2=[sorted(arr3)[x] for x in arr1]       # det(A)= (-1)^n * Sum[(-1)^alpha * f(h)]
    count=0
    for x in range(n):
        if arr3[x]!=arr2[x]:
            (arr2[arr2.index(arr3[x])],arr2[x])=(arr2[x],arr3[x])
            count+=1
    return count

multiplier = 1 if n%2 else -1
for x,sub in enumerate(subs_S):
    det=0
    for c in sub:
        p=no_of_ways(sorted_rows,c)            # alpha in the formula for det(A)
        sign=(-1)**(p)
        for posn in list(zip(sorted_rows,c)):
           sign*=VA[posn]            
        det+=sign
    if x!=n:det=((-1)**x)*multiplier*det    
    Det_subs.append(det)
   
det_coeff=Det_subs[n]                          # Determinant of coefficient matrix
for i in range(n):
    print("x%d ="%(i+1),Det_subs[i]/det_coeff,end=", " if i!=(n-1) else ' ')
    
    
