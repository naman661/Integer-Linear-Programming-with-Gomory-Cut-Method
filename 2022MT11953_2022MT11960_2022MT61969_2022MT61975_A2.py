import copy
import math
import numpy as np
from numpy import linalg as LA # import linear algebra 
from numpy.linalg import inv

def row_operations(A,pivot_row_index,pivot_column_index):

    # pivot element = k
    k = A[pivot_row_index][pivot_column_index]

    # new matrix B
    B = copy.deepcopy(A)

    n = len(A[0])
    m = len(A)

    for i in range(n):
        B[pivot_row_index][i] = B[pivot_row_index][i]/k

    for i in range(m):

        if i == pivot_row_index:
            continue

        multiplication_factor = A[i][pivot_column_index]
        for j in range(n):
            B[i][j] = B[i][j] - multiplication_factor*B[pivot_row_index][j]

    for i in range(len(B)):
       for j in range(len(B[0])):
          B[i][j] = round(B[i][j],14)

    return B

def iteration(A1,basis_indices):

    A = copy.deepcopy(A1)
    for i in range(len(A)):
       for j in range(len(A[0])):
          A[i][j] = round(A[i][j],14)


    m = len(A) 
    n = len(A[0]) 
    dual_unbounded = False
    optimal = False

    while(not dual_unbounded and not optimal):

        #  checking if optimal
        tempo = True
        for i in range(1,m):
            if A[i][0]<0:
                tempo = False
                break
        if tempo:
            optimal = True

        if optimal:
            continue



        ## if not optimal
        
        # searching for pivot row (smallest index row)
        pivot_row_index = 0
        for i in range(1,m):
            if(A[i][0]<0):
                pivot_row_index = i
                break
        
        # checking if dual unbounded
        tempo1 = True
        for i in range(1,len(A[0])):
            if A[pivot_row_index][i]<0:
                tempo1 = False
                break
        
        if tempo1:
            dual_unbounded = True
        
        if dual_unbounded:
            continue



        ## if not unbounded
        
        # searching for pivot column (using Blands rule)
        current_min = 10000000000
        pivot_cloumn_index = -1
        for i in range(1,len(A[0])):
            if A[pivot_row_index][i] < 0:
                tempo = A[0][i]/abs(A[pivot_row_index][i])
                if tempo<current_min:
                    current_min = tempo
                    pivot_cloumn_index = i

        # updating basis_indices
        basis_indices[pivot_row_index-1] = pivot_cloumn_index - 1


        # applying row transformations
        A = row_operations(A,pivot_row_index,pivot_cloumn_index)
    
    if(dual_unbounded):
        return "dual_unbounded",[]

    else:
        return A,basis_indices

def toStandard(A,b,c,constraint_type):
    A2 = []
    b2 = []
    c2 = []

    s=[]
    ns = []

    for i in range(len(A)):
        if constraint_type[i] != 0:
            s.append(i)
        else:
            ns.append(i)
        
  
    for i in s:
       A2.append(A[i])
    for i in ns:
       A2.append(A[i])

    for i in range (len(s)):
       b2.append([b[s[i]][0]])
       for j in range(len(s)):
          if i == j:
             if constraint_type[s[i]] == -1:
                A2[i].append(-1)
             else:
                A2[i].append(1)
          else:
             A2[i].append(0)
    
    for i in range(len(s),len(A)):
       for j in range(len(s)):
          A2[i].append(0)

    for i in ns:
       b2.append([b[i][0]])

    c2 = c.copy()
    for i in range(len(s)):
       c2.append(0)
    
    for i in range(len(A2)):
       if (b2[i][0]<0):
        b2[i][0] *= -1
        A2[i] = [x * -1 for x in A2[i]]

    # row = len(A2)

    # for i in range(row):
    #    c2.append(0)
    #    for j in range(row):
    #       if i==j:
    #          A2[i].append(1)
    #       else:
    #          A2[i].append(0)
    return A2,b2,c2 
  
def simplex_iteration(A, b, C, m: int, n: int,is_auxillary,Answer_dictionary,Global_B,Global_Basis):
    if(is_auxillary==-69):
      return
    #intialization
    Iteration=0
    Z=0
    X=np.zeros((n+m))
    XB=np.zeros((m))
    CB=np.zeros((m))
    XN=np.zeros((n))
    CN=np.zeros((n))
    RC = np.zeros((n+m))
    Basis:int=np.zeros((m))
    B = np.zeros((m,m))
    NB = np.zeros((m,n))
    Index_Enter=-1
    Index_Leave=-1
    eps = 1e-12

    for i in range(0,m):
        Basis[i]=n+i
        for j in range(0,m):
         B[i, j]=A[i,n+j]
        for j in range(0,n):
         NB[i, j]=A[i,j]

    for i in range(0,n):
        CN[i]=C[i]
    if(is_auxillary==0):
      Basis = Global_Basis
      B = Global_B
    for i in range(len(Basis)):
      CB[i] = C[int(Basis[i])]
    RC=C-np.dot(CB.transpose(),np.dot(inv(B),A))
    MaxRC=0
    X=np.dot(inv(B),b)
    Z=np.dot(CB,X)
    for i in range(0,n+m):
        if(MaxRC<RC[i]):
         MaxRC=RC[i]
         Index_Enter=i
    Down_Table_initial = np.concatenate((b,A),1)
    Down_Table_initial = np.dot(inv(B),Down_Table_initial)

    if(is_auxillary!=1):  
      Answer_dictionary["initial_tableau"]=(np.round(Down_Table_initial,14))

    while(MaxRC > eps):
      Iteration=Iteration+1

      Index_Leave=-1
      MinVal=1000000

      for i in range(0,m):
       if(np.dot(inv(B),A)[i,  Index_Enter] > 0):
         bratio=np.dot(inv(B),b)[i]/np.dot(inv(B),A)[i,  Index_Enter]

         if(MinVal > bratio ):
           Index_Leave=i

           MinVal=bratio


      if (Index_Leave == -1):

        Down_Table_final = np.concatenate((b,A),1)
        Down_Table_final = np.dot(inv(B),Down_Table_final)

        Answer_dictionary["final_tableau"]=(np.round(Down_Table_final,14))

        Answer_dictionary["solution_status"]="unbounded"
        return Z,X,RC
      
      Basis[Index_Leave]=Index_Enter 

      for i in range(m-1,0,-1):
        if(Basis[i] < Basis[i-1]):
            temp=Basis[i-1]
            Basis[i-1]=Basis[i]
            Basis[i]=temp

      for i in range(0,m):
          for j in range(0,n+m):
              if(j==Basis[i]):
                B[:, i]=A[:,j]
                CB[i]=C[j]

      RC=C-np.dot(CB.transpose(),np.dot(inv(B),A))
      MaxRC=0
      for i in range(0,n+m):
        if(MaxRC<RC[i]):
         MaxRC=RC[i]
         Index_Enter=i

      X=np.dot(inv(B),b)
      Z=np.dot(CB,X)

    Down_Table_final = np.concatenate((b,A),1)
    Down_Table_final = np.dot(inv(B),Down_Table_final)

    if(is_auxillary==0):     

      Answer_dictionary["final_tableau"]=(np.round(Down_Table_final,14))

      Answer_dictionary["solution_status"]="optimal"

      ans = [0]*len(RC)
      for i in range(len(Basis)):
        ans[int(Basis[i])] = np.round(X[i][0],14)

      Answer_dictionary["optimal_solution"]=ans

      Answer_dictionary["optimal_value"]=np.round(Z,14)[0]
      Answer_dictionary["Final_B"]=B
      Answer_dictionary["Final_Basis"]=Basis
    else:
      if(Z<0):
        Answer_dictionary["solution_status"]="infeasible"
        return B,Basis,-69
      else:
        return B,Basis,0

def simplex_algo():
  Answer_dictionary = {"initial_tableau":np.array([[]]),"final_tableau":np.array([[]]),"solution_status":"","optimal_solution":[],"optimal_value":0,"Final_B":np.array([[]]),"Final_Basis":np.array([]),"standard_A":np.array([[]]),"standard_B":np.array([[]]),"standard_C":np.array([]),"is_max":True,"no_of_var":0}
  is_max = True
  A1 = []
  b1 = []
  constraint_type = []
  C1 = []
  file = open("input_ilp.txt", "r")
  content=file.readlines()

  file.close()
  if(content[1][:-1]!="maximize"):
    is_max = False

  index_A = -1
  index_b = -1
  index_C = -1
  index_constraint = -1
  for i in range(len(content)):
    if(content[i][:-1]=="[A]"):
      index_A = i
    elif(content[i][:-1]=="[b]"):
      index_b = i
    elif(content[i][:-1]=="[constraint_types]"):
      index_constraint = i
    elif(content[i][:-1]=="[c]"):
      index_C= i

  for i in range(index_A+1,index_b-1): 
    temp = [int(x) for x in content[i][:-1].split(',')]
    A1.append(temp)

  for i in range(index_b+1,index_constraint-1): 
    temp = [int(x) for x in content[i][:-1].split(',')]
    b1.append(temp)  

  for i in range(index_constraint+1,index_C-1): 
    if(content[i][:-1]=="<="):
      constraint_type.append(1)
    elif(content[i][:-1]==">="):
      constraint_type.append(-1)
    else:
      constraint_type.append(0) 

  C1 = [int(x) for x in content[index_C+1][:].split(',')] 

  if(is_max==False):
    C1 = [x * -1 for x in C1]
    Answer_dictionary["is_max"] = False
  Answer_dictionary["no_of_var"] = len(A1[0])
  A1,b1,C1 = toStandard(A1.copy(),b1.copy(),C1.copy(),constraint_type)  
  
  temp_A1 = copy.deepcopy(A1)
  temp_b1 = copy.deepcopy(b1)
  temp_C1 = copy.deepcopy(C1)

  temp_A1 = np.array(temp_A1)
  temp_b1 = np.array(temp_b1)
  temp_C1 = np.array(temp_C1)

  Answer_dictionary["standard_A"] = temp_A1
  Answer_dictionary["standard_B"] = temp_b1
  Answer_dictionary["standard_C"] = temp_C1
  

  var_to_denote_feasible = 1
  C_temp = [0]*len(A1[0])+[-1]*len(constraint_type)
  A=np.array(A1.copy())
  b=np.array(b1.copy())
  C1 = np.array(C1)
  C3=C_temp.copy()
  C3 = np.array(C3)
  A3=A1.copy()
  for i in range(len(constraint_type)):
    for j in range((len(constraint_type))):
        if(i==j):
          A3[i].append(1)
        else:
          A3[i].append(0)
  A3 = np.array(A3)
  b3=np.array(b1.copy())
  Global_Basis = np.array([])
  Global_B = np.array([[]])
  Global_B,Global_Basis,var_to_denote_feasible = simplex_iteration(A3,b3,C3,len(A3),len(A3[0])-len(A3),var_to_denote_feasible,Answer_dictionary,Global_B,Global_Basis)
  simplex_iteration(A,b1,C1,len(A),len(A[0])-len(A),var_to_denote_feasible,Answer_dictionary,Global_B,Global_Basis)
  if(is_max==False):
    Answer_dictionary["optimal_value"] = Answer_dictionary["optimal_value"]*-1
  return Answer_dictionary

def cutting_plane(A1,basis_indices):

    ans_dict = {"solution_status":"","final_solution":[],"number_of_cuts":0,"optimal_value":0}
    
    A = copy.deepcopy(A1)

    m1 = len(A)
    n1 = len(A[0])

    optimal = False
    dual_unbounded = False
    
    eps = 1e-4
    # assuming primal is bounded

    count = 0
    n = len(A[0])
    while(not optimal and not dual_unbounded):
        
        count += 1
        print(count)
        if count == 100:
           optimal = True
           break

        m = len(A)
        n = len(A[0])

        # rounding off
        for i in range(len(A)):
            for j in range(len(A[0])):
                if(abs(A[i][j]-round(A[i][j]))<eps):
                    A[i][j] = round(A[i][j])

        tempo = True
        for i in range(1,m):
            tempo1 = A[i][0] + 0.0
            if not tempo1.is_integer():
                tempo = False
                break
        if tempo:
            optimal = True
    
        if optimal:
            continue

        
        ## if not optimal

        # searching non integer A[i][0] (source row)

        source_row_index = 0

        for i in range(m):
            tempo1 = A[i][0] + 0.0
            if not tempo1.is_integer():
                source_row_index = i
                break
        
        
        # creating row to add
        r = [0.0]*(n-1)

        for j in range(n-1):
            if j not in basis_indices:
                r[j] = -1 * ((A[source_row_index][j+1] + 0.0)%1)

        r.append(1)
        r.insert(0,-1 * ((A[source_row_index][0]+0.0)%1))

        # creating column in A
        for i in range(m):
            A[i].append(0)

        # adding row to A
        A.append(r)

        # adding new index to basis_indices
        basis_indices.append(len(A[0]) - 2)


        
        
        # performing dual-simplex iteration
        print(count)
        A,basis_indices = iteration(A,basis_indices)

        print(count)
        #checking if dual_unbounded
        if A == "dual_unbounded":
            dual_unbounded = True
        else:
            for i in range(len(A)):
                for j in range(len(A[0])):
                    if(abs(A[i][j]-round(A[i][j]))<eps):
                        A[i][j] = round(A[i][j])
            for i in range(len(A)):
                for j in range(len(A[0])):
                    A[i][j] = round(A[i][j],14)
        
    if dual_unbounded:
        ans_dict["solution_status"] = "infeasible"
        ans_dict["number_of_cuts"] = n-n1
        return ans_dict

    if optimal:
        ans_dict["solution_status"] = "optimal"

        m = len(A)

        sol = [0.0]*(n-1)
        for i in range(1,m):
            sol[basis_indices[i-1]] = A[i][0]

        ans = sol[:n1-1]

        for i in range(len(ans)):
           ans[i] = round(ans[i])
        ans_dict["final_solution"] = ans

        ans_dict["optimal_value"] = round(-1*A[0][0])

        ans_dict["number_of_cuts"] = n-n1

    return ans_dict

def helper():

    D = simplex_algo()
    no_of_variables = D["no_of_var"]

    # checking if primal is unbounded
    if D["solution_status"] == "unbounded":
        print("initial_solution: ")

        print("final_solution: ")

        print("solution_status: unbounded")

        print("number_of_cuts: ")

        print("optimal_value: ")
        return
    
    # checking if primal is infeasible
    if D["solution_status"] == "infeasible":
        print("initial_solution: ",)

        print("final_solution: ")

        print("solution_status: infeasible")

        print("number_of_cuts: ")

        print("optimal_value: ")
        return
    

    # primal is optimal
    A1 = D["standard_A"]
    b1 = D["standard_B"]
    c1 = D["standard_C"]
    c1 = -1*c1


    basis_indices = D["Final_Basis"]
    basis_indices = [ int(x) for x in basis_indices ]
    basis = D["Final_B"]
    initial_solution = D["optimal_solution"]
    initial_solution_value = D["optimal_value"]
    is_max = D["is_max"]
    if is_max:
       initial_solution_value = -1 * initial_solution_value

    dimensions = A1.shape
    m,n = dimensions

    inv_basis = np.linalg.inv(basis)
    inv_B_mul_A = np.dot(inv_basis,A1)


    CB = [[]]
    for i in range(m):
        CB[0].append(c1[basis_indices[i]])

    CB = np.array(CB)


    #converting c1 to matrix form
    C = [[]]
    for i in range(n):
        C[0].append(c1[i])
    C = np.array(C)
    c1 = C
    

    # reduced costs
    RC = c1 - np.dot(CB,inv_B_mul_A)
    #creating first row
    r1 = []
    r1.append(-1*initial_solution_value)
    for i in range(n):
        r1.append(RC[0][i])
    
    # creating initial tableau
    
    Initial_tableau = []
    Initial_tableau.append(r1)

    for i in range(m):
        r = []
        r.append(initial_solution[basis_indices[i]])
        for j in range(n):
            r.append(inv_B_mul_A[i][j])
        Initial_tableau.append(r)

    for i in range(len(Initial_tableau)):
       for j in range(len(Initial_tableau[0])):
          Initial_tableau[i][j] = round(Initial_tableau[i][j],14)


    

    ans_dict = cutting_plane(Initial_tableau,basis_indices)

    solution_status = ans_dict["solution_status"]

    if solution_status == "infeasible":

        print("initial_solution: ",end="")
        for i in range(no_of_variables-1):
            print(str(initial_solution[i]) + ",",end=" ")
        print(initial_solution[no_of_variables-1])

        print("final_solution: ")

        print("solution_status: infeasible")

        number_of_cuts = ans_dict["number_of_cuts"]
        print("number_of_cuts: " + str(number_of_cuts))

        print("optimal_value: ")
    
    if solution_status == "optimal":
       
        print("initial_solution: ",end="")
        for i in range(no_of_variables-1):
            print(str(initial_solution[i]) + ",",end=" ")
        print(initial_solution[no_of_variables-1])


        final_solution = ans_dict["final_solution"]
        print("final_solution: ",end="")
        for i in range(no_of_variables-1):
            print(str(final_solution[i]) + ",",end=" ")
        print(final_solution[no_of_variables-1])


        print("solution_status: optimal")


        number_of_cuts = ans_dict["number_of_cuts"]
        print("number_of_cuts: " + str(number_of_cuts))


        optimal_value = ans_dict["optimal_value"]
        if is_max:
           optimal_value = -1 * optimal_value
        print("optimal_value: " + str(optimal_value))

def gomory_cut_algo():
  helper()

gomory_cut_algo()

