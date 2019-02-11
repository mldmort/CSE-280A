import numpy as np
from methods import Node
from sets import Set

if __name__ == "__main__":

  n = 10
  tetha = 1e-5
  mu = 1e-8
  l = int(1e6)
  N = int(1e4)

  Node_list = [[]]
  for i in xrange(n):
    Node_list[0].append(Node(i+1))

  for i in xrange(n,1,-1):
    n_list = range(i)
    Node_list.append([])
    rand2 = np.random.choice(n_list,2,replace=False)
    Node_list[-1].append(Node())
    Node_list[-1][-1].left = Node_list[-2][rand2[0]]
    Node_list[-1][-1].right = Node_list[-2][rand2[1]]
    Node_list[-1][-1].nch = 2

    n_list.remove(rand2[0])
    n_list.remove(rand2[1])

    for i in n_list:
      Node_list[-1].append(Node())
      Node_list[-1][-1].left = Node_list[-2][i]
      Node_list[-1][-1].nch = 1
    
  #for i in range(len(Node_list)):
  #  for j, node in enumerate(Node_list[i]):
  #    dec = []
  #    dec = node.getDec(dec)
  #    print '> level: ' , i , ' ,node #: ' , j , ' ,dec: ' , dec
  #my_dec = Node_list[4][0].getDec()
  
  alpha = 0.0
  #mutation_list = range(l)
  mutation_set = Set([])
  #for k in xrange(2,n+1):
  for k in xrange(2,2+1):
    Tk = int(np.exp(-alpha*(n-k+1))*N*2/k/(k-1))
    avg_mutation = Tk * mu * l
    n_mutation = np.random.poisson(avg_mutation, k)
    for ik in xrange(k):
      my_mutation = 0
      while my_mutation < n_mutation[ik]:
        myrand = np.random.randint(0,l-1)
        if (not(myrand in mutation_set)):
          mutation_set.add(myrand)
          my_mutation += 1

    #for kk in xrange(k):
    #  mutation_k = np.random.choice(mutation_list,n_mutation[kk],replace=False)
    #  for im in mutation_k:
    #    mutation_list.remove(im)
    
    
    print 'Tk: ' , Tk , ' avg_mutation: ' , avg_mutation , ' n_mutation: ' , n_mutation
    #print 'mutation_k: ' , mutation_k
  
  #print 'rand2: ' , rand2
  #print 'n_list: ', n_list
  #print 'Node_list: ', Node_list   
  print 'Node_list size: ', len(Node_list) 
  #for i in xrange(n): 
  #  print len(Node_list[i])

  #print 'my_dec: ' , my_dec
