import numpy as np
from numpy import linalg as LA # import linear algebra 
from numpy.linalg import inv
import sys

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
        # print("CN: ", CN[i]) 
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

    # print("Basis", Basis)
    # C_initial = C.copy()
    # C_initial = np.insert(C_initial,0,0)
    Down_Table_initial = np.concatenate((b,A),1)
    Down_Table_initial = np.dot(inv(B),Down_Table_initial)
    # Initial_Tableau = np.concatenate((C_initial.reshape(1,-1),Down_Table_initial),0)
    if(is_auxillary!=1):  
      # print("----------------Initial_Tableau----------------")
      # print(np.round(Initial_Tableau,2))
      Answer_dictionary["initial_tableau"]=(np.round(Down_Table_initial,2))
      # print("-----------------------------------------------")

    while(MaxRC > eps):
      Iteration=Iteration+1
      # print("=> Iteration: ",Iteration)

      # print(" Index_Enter: ",  Index_Enter)
      Index_Leave=-1
      MinVal=1000000
      # print("Enter B: ",B)
      for i in range(0,m):
       if(np.dot(inv(B),A)[i,  Index_Enter] > 0):
         bratio=np.dot(inv(B),b)[i]/np.dot(inv(B),A)[i,  Index_Enter]
        #  print("  bratio: ", bratio)
         if(MinVal > bratio ):
           Index_Leave=i
          #  print("  Index_Leave: ",Index_Leave)
           MinVal=bratio
          #  print("  MinVal: ", MinVal)

      if (Index_Leave == -1):
        # C_final = RC.copy()
        # C_final = np.insert(C_final,0,-Z[0])
        Down_Table_final = np.concatenate((b,A),1)
        Down_Table_final = np.dot(inv(B),Down_Table_final)
        # Final_Tableau = np.concatenate((C_final.reshape(1,-1),Down_Table_final),0)
        # print("----------------Final_Tableau----------------")
        # print(np.round(Final_Tableau,2))
        Answer_dictionary["final_tableau"]=(np.round(Down_Table_final,2))
        # print("-----------------------------------------------")
        Answer_dictionary["solution_status"]="unbounded"
        return Z,X,RC
      
      Basis[Index_Leave]=Index_Enter 
      # print("before updated Basis", Basis)
      # print("  Index_Leave: ",Index_Leave)
      for i in range(m-1,0,-1):
        if(Basis[i] < Basis[i-1]):
            temp=Basis[i-1]
            Basis[i-1]=Basis[i]
            Basis[i]=temp

      # print("updated Basis", Basis)

      for i in range(0,m):
          for j in range(0,n+m):
              if(j==Basis[i]):
                B[:, i]=A[:,j]
                CB[i]=C[j]

      # print("Exit Basis", Basis)
      # print("Exit B: ",B)

      RC=C-np.dot(CB.transpose(),np.dot(inv(B),A))
      MaxRC=0
      for i in range(0,n+m):
        if(MaxRC<RC[i]):
         MaxRC=RC[i]
         Index_Enter=i
      # print("MaxRC",MaxRC)
      X=np.dot(inv(B),b)
      Z=np.dot(CB,X)
    # C_final = RC.copy()
    # C_final = np.insert(C_final,0,-1*Z)
    Down_Table_final = np.concatenate((b,A),1)
    Down_Table_final = np.dot(inv(B),Down_Table_final)
    # Final_Tableau = np.concatenate((C_final.reshape(1,-1),Down_Table_final),0)
    if(is_auxillary==0):     
      # print("----------------Final_Tableau----------------")
      # print(np.round(Final_Tableau,2))
      Answer_dictionary["final_tableau"]=(np.round(Down_Table_final,2))
      # print("-----------------------------------------------")
      Answer_dictionary["solution_status"]="optimal"
      # print("optimal")
      ans = [0]*len(RC)
      for i in range(len(Basis)):
        if (Basis[i]<n):
          ans[int(Basis[i])] = np.round(X[i][0],2)
      # print("-----------------------------------------------")
      Answer_dictionary["optimal_solution"]=ans
      # print("Optimal Solution: ",ans)
      # print("-----------------------------------------------")
      # print("Optimal Value: ",np.round(Z,2)[0])
      Answer_dictionary["optimal_value"]=np.round(Z,2)[0]
      Answer_dictionary["Final_B"]=B
      Answer_dictionary["Final_Basis"]=Basis
    else:
      if(Z>0):
        Answer_dictionary["solution_status"]="infeasible"
        return -69
      else:
        return B,Basis,0

def simplex_algo():
  Answer_dictionary = {"initial_tableau":np.array([[]]),"final_tableau":np.array([[]]),"solution_status":"","optimal_solution":[],"optimal_value":0,"Final_B":np.array([[]]),"Final_Basis":np.array([]),"standard_A":np.array([[]]),"standard_B":np.array([[]]),"standard_C":np.array([]),"is_max":True}
  is_max = True
  A1 = []
  b1 = []
  constraint_type = []
  C1 = []
  file = open("input.txt", "r")
  content=file.readlines()
  # print(content)
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

  A1,b1,C1 = toStandard(A1.copy(),b1.copy(),C1.copy(),constraint_type)  

  Answer_dictionary["standard_A"] = A1
  Answer_dictionary["standard_B"] = b1
  Answer_dictionary["standard_C"] = C1

  print("A1",A1)
  print("b1",b1)
  print("C1",C1)

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

ANSW = simplex_algo()
print(ANSW)