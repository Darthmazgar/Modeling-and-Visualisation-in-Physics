import time
P=input()
N=range(20)
while 1:
 for i in N:print''.join(' *'[i*20+j in P]for j in N)
 time.sleep(.1);Q=[(p+d)%400 for d in(-21,-20,-19,-1,1,19,20,21)for p in P];P=set(p for p in Q if 2-(p in P)<Q.count(p)<4)
