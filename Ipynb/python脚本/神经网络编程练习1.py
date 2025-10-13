#!/usr/bin/env python
# coding: utf-8

# In[1]:


for n in range(10):
    print(n)
    pass
print("done")


# In[2]:


# the following prints out the cube of 2
print(2**3)


# In[3]:


import numpy


# In[4]:


a = numpy.zeros( [3,2] )
print(a)


# In[5]:


import matplotlib.pyplot
get_ipython().run_line_magic('matplotlib', 'inline')


# In[6]:


a[0,0] = 1
a[0,1] = 2
a[1,0] = 9
a[2,1] = 12
print(a)


# In[7]:


matplotlib.pyplot.imshow(a, interpolation="nearest")


# In[8]:


for n in range(100):
    print("The square of",n,"is",n*n)
    pass
print("done")


# In[ ]:




