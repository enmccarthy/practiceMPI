import numpy as np
import h5py

f = h5py.File('/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5/21688988/univ_ics_2019-03_a9830598.hdf5', 'r')
print(list(f.keys()))

#x = np.load('/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/dim512_cube_nt4_npz/train/train_a1114020_int16.npz')
# /p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/dim512_cube_nt4_npz/train/train_a1114020_int16.npz
data = f[u'full']
dataP = np.array(data)
print(data.shape)
print(data[:256][:256][0][0])
dataR = dataP.reshape(4,512,512,512)
pt1 = (dataP[:256,:256,:,:])
print(pt1.shape)
file1 = open('output0.txt', 'r')
alldata = file1.read()
arr = alldata.split(" ");
for ind, num in enumerate(np.ndarray.flatten(pt1)):
    currRead = int(arr[ind])
    #print(currRead, num)
    if not currRead == num: 
        print(ind)
        print(currRead)
        print(num)
        break
print("DONE")
