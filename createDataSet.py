import numpy as np
import h5py

d1 = np.random.random_integers(0,100,size=(4,512,512,512))
print(d1[0].shape)
print(d1.dtype)
d1 = d1.astype(np.int16)
hf = h5py.File('/p/gpfs1/emccarth/test10.h5', 'w')
hf.create_dataset('dataset_1', data=d1)
hf.close()
