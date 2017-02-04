# Computing determinant and permanent of a square sparse matrix by Numerical Structure approach
# Assumptions: The square matrix is non-singular (i.e. Determinant exists)

m=int(input("Enter the order 'm' of your square matrix: "))
print("Input the %d x %d matrix: "%(m,m))
elements=dict()     # All non zero entries

# Storing only the non-zero entries in the form "elements[(row,column)] = value" as a dictionary
for i in range(m):
    for (j,v) in enumerate(list(map(int,input().split()))):
        if v!=0:elements[(i,j)]=v

# Numerical structure 'S' for the matrix
S=dict()
for row,col in elements:S[row]=[]       # Initializing

for row,col in sorted(elements):
    S[row].append(col)

sorted_S=list()                         # The rows of S are sorted in increasing length of set H
sorted_rows=list()                      # row indices after sorting

for row_ind,H in sorted(enumerate(S.values()),key=lambda x:len(x[1])):   
    sorted_S.append(H)                  # H -> set of column indices of non zero entries in each row in S
    sorted_rows.append(row_ind)

# Taking permutations of each H such that column indices are unique    
unique_cols=[]
for i in sorted_S[0]:unique_cols.append([i])

for h in range(1,m):
    temp=len(unique_cols)
    for y in range(len(unique_cols)):
        for x in sorted_S[h]:
            if x not in unique_cols[y]:
                unique_cols.append(unique_cols[y]+[x])
    unique_cols=list(unique_cols[temp:])

# Computing Determinant and Permanent
det=0					                         # Stores Determinant
perm=0					                       # Stores Permanent

def no_of_ways(arr1,arr2):        		 # Interchanges needed
    arr3=arr2[:]
    count=0
    for x in range(m):
        if arr1[x]!=arr3[x]:
            arr3[arr3.index(arr1[x])],arr3[x] = (arr3[x],arr1[x])
            count+=1
    return count

for c in unique_cols:
    p=(m - no_of_ways(sorted_rows,c))
    sig=(-1)**(p)
    const=1
    for posn in list(zip(sorted_rows,c)):
        sig*=elements[posn]
        const*=elements[posn]    
    perm+=const
    det+=sig
det=((-1)**m)*det

print("Determinant of the matrix by numerical structure approach is",det)
print("Permanent of the matrix by numerical structure approach is",perm)
        
                
                                         
    

    
