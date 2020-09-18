


from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


# In[26]:


data = pd.read_csv('cmake-build-debug/predict_data.csv', header=None)
learn_data = pd.read_csv('cmake-build-debug/learn_data.csv', header=None)


# In[27]:

plt.plot(data[0] + 1872, data[1])
plt.plot(learn_data[0] + 1872, learn_data[1])
plt.show()


# In[28]:


# ## Part 2

# In[29]:


data = pd.read_csv('cmake-build-debug/predict_data2.csv', header=None)
learn_data = pd.read_csv('cmake-build-debug/learn_data2.csv', header=None)


# In[30]:


plt.plot(data[0] + 1872, data[1])
plt.plot(learn_data[0] + 1872, learn_data[1])
plt.show()

# In[31]:


# In[ ]:




