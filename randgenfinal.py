import random
import time
import numpy as np

el1=raw_input("elemento 1 e numero 1 ")
n1=input()
n1=int(n1)
l=[el1]*n1
el2=raw_input("elemento 2 e numero 2 ")
n2=input()
n2=int(n2)
l2=[el2]*n2
for i in range(len(l2)):
    l.append(l2[i])
#print(l1)

print(l)

with open('Cu55','a+') as f:   #change file here; index guide in line 31  
    lines = f.readlines()
f.close()
lines = [w.replace('\t', '') for w in lines]
lines = [w.replace('\n', '') for w in lines]
y=0


start = time.time()
for i in range(1000):
        lr = random.sample(l, len(l))
        l=lr
        y=y+1    #change 1=1,2=57,3=113,4=169,5=225,6=281,7=337,8=393,9=449,10=505... if using several generating geometries (one at a time)
        yy=str(y)
        new= open(str(el1)+str(n1)+str(el2)+str(n2)+'_'+yy+'.xyz', 'a+')     ##DEFINIR O NOME DOS ARQUIVOS A SEREM GERADOS
        k=np.array(l)
        #k=l        
        new.write(str(55)+'\n'+'\n') #write 8, skip line twice
        for i, r in zip(k,lines):
                new.write(i+'   '+r+'\n') #write element + coordinates
        f.close()

#        print(str(yy))

end = time.time()
tot=end - start
#print(str(tot))
#print(lr)

