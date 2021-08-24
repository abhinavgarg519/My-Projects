#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


# In[2]:


L= (10)**(-6)#lenth of cell
a= 2.5*(10**-10)#burger vector or spacing between points
D= int(L/a)+1 #No. of positions in the cell
n=50 #No. of dislocations
G=5.5*(10**10)
v=0.324
density = n/L**2
stress=0.3*(10**9) #Externally applied shear stress sigma-xy 
print(L,a,D,n,G,v,density)


# In[3]:


x=np.arange(0,D) 
print(x)
b=np.zeros(D)
b[0:int(D/2)]=1
b[int(D/2):]=-1
import random
random.shuffle(b)
print(b)


# In[4]:


import random
y=np.zeros(D)
# Putting n random values between[1,D] in y representing n dislocations
i= random.sample(list(np.arange(1,D+1)),n)
j= random.sample(list(x),n)
for p in range (0,n):
    k=j[p]
    y[k]=i[p]
print(y)
np.unique(y).shape
#Don't Panic, it includes zero also.


# In[5]:


initloc= pd.DataFrame((x,y,b),['X','Y','B'])
#Making a data frame containing positions of dislocations and the sign of their burger vector 
#Here they are asigned to a index which represents positions for dislocations in the structure.
initloc= initloc.T
initloc.loc[:,'B'][initloc['Y']==0]=0 #Assigning "burger vectors sign" to zero for rows in which there is no dislocation.
print(initloc[initloc['Y']!=0].count())
#count of rows where y is not equal to zero for confirming number of dislocations
print(initloc)


# In[6]:


plt.scatter(initloc[initloc['B']>0]['X'],initloc[initloc['B']>0]['Y'],marker='o',c='b',linewidths=.1)
plt.scatter(initloc[initloc['B']<0]['X'],initloc[initloc['B']<0]['Y'],marker='o',c='r',linewidths=.1)
plt.xlim(0,4000)
plt.ylim(0,4000)
plt.savefig("Initial Location.png")


# In[7]:


cond = initloc[initloc['B']!=0] #separate dataframe in which there are only dislocations and not all the points
ind = cond.index #index of cond dataframe act as unique identity of dislocation and will be used to iterate through each dislocation
print(ind)
print(cond.loc[ind[0],'X']) #example of how index can be used 
cond.to_csv("Initital Location.csv") #saving inital location to a csv file
print(cond)


# In[8]:


count=0
fx= pd.DataFrame()#storing force data
pos=pd.DataFrame()#storing position data
pos =pos.append(cond)
print(pos)


# In[54]:


#force code
c=G/(2*math.pi*(1-v))
dxmax = 200*a
ntstep= 200#no. of time steps

for m in range (count,(count+ntstep)):#instead of taking(0,ntstep) I have taken this so that for each timestep has a unique no.
    n=len(cond.index)
    fmax=0
   
    for i in  range (0,n):
        fx.loc[i,m]=0
        for j in range(0,n):
            if i!=j:
                dx= (cond.iloc[i]['X']-cond.iloc[j]['X'])*a #x[i]-x[j]
                dy= (cond.iloc[i]['Y']-cond.iloc[j]['Y'])*a #y[i]-y[j]
                dsq= (dx**2)+(dy**2)
                dist= dsq**0.5
                if cond.iloc[i]['Y']==cond.iloc[j]['Y'] and dist<=(6*a):#Dislocation anihilation
                    cond.loc[ind[i],'B']=0
                    cond.loc[ind[j],'B']=0
                    fx.loc[i,m]=0
                    fx.loc[j,m]=0
                
                ffx= (cond.iloc[i]['B'])*(cond.iloc[j]['B'])*a*a*c*dx*(dx**2-dy**2)/(dsq**2)#calculating force
                fx.loc[i,m]= fx.loc[i,m]+ffx
        fx.loc[i,m]=fx.loc[i,m]+(stress*a*(cond.iloc[i]['B']))#force due to external shear stress

        if abs(fx.loc[i,m])>fmax:
            fmax= abs(fx.loc[i,m])
    tstep = dxmax/fmax#size of time step
    
    for q in range(0,n):
        row=ind[q] #for each q ind[q] will give us unique no. of each dislocation which we can use to store and edit info.
        cond.loc[row,'X']=cond.loc[row,'X']+int(fx.loc[q,m]*tstep/a) #updating the x-coordinates of dislocations in integers
        #used int to remove movement across periodic boundaries
        if cond.loc[row,'X']>(D-1):
            cond.loc[row,['X']]-=D
        if cond.loc[row,'X']<0:
            cond.loc[row,['X']]+=D
        pos.loc[row,m]=0
        pos.loc[row,m]=cond.loc[row,'X']
    
    cond= cond[cond['B']!=0]#Removed the dislocation which are anihilated
    ind=cond.index
count+=ntstep


# In[55]:


print(fmax, tstep)
print(cond)
cond.to_csv("final location.csv")
print(cond.count())#counting the no. of dislocations


# In[56]:


print(cond.count())


# In[57]:


print(fx.shape)


# In[58]:


plt.scatter(cond[cond['B']>0]['X'],cond[cond['B']>0]['Y'],marker='o',c='b',linewidths=.1)
plt.scatter(cond[cond['B']<0]['X'],cond[cond['B']<0]['Y'],marker='o',c='r',linewidths=.1)
plt.xlim(0,4000)
plt.ylim(0,4000)
plt.savefig("Final location.png")


# In[59]:


print(count)


# In[60]:


fx.to_csv("funcvalue.csv")#contains force value for each dislocation(rows) for each timestep(columns)


# In[61]:


pos.to_csv("posvalue.csv")#contains position of dislocations from initial to final timestep


# In[62]:


print(pos.shape)


# In[ ]:





# In[ ]:




