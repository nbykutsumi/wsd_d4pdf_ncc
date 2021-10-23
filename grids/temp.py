# %%
import numpy as np
from numpy import ma
#%matplotlib inline
from myfunc.grids import expand_map_global_2d, shift_map_global, karnel_pooling_map2D_global, shift_map_regional,karnel_pooling_map2D_regional


a=np.arange(1,10).reshape(3,3)


print(a)
print()
#b=shift_map_regional(a, -1,0)
b= karnel_pooling_map2D_regional(a,1,1,func="sum")
print("---------------")
print(b)

