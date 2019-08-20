import numpy as np
x = np.load('/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/dim512_cube_nt4_npz/train/train_a1114020_int16.npz')
# /p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/dim512_cube_nt4_npz/train/train_a1114020_int16.npz
data = x['data']
print(data[0][0][0][:256])
print(data[0][0][0][0][256:])
pt1 = (data[0,:,:,256:,:256])
print(pt1.shape)
file1 = open('output2.txt', 'r')
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
