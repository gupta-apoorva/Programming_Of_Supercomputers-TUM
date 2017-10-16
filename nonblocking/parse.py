import re
import numpy as np

Size = [64,512,1024,2048,4096,8192]
p=[8,16,32,64]
Pros = [8,24,56,120]

plot_ind=4

num_file_to_merge = 10
data=np.zeros((num_file_to_merge, 6, 120, 5))
for i in range(num_file_to_merge):
    hand = open("out_origin"+str(i+1)+".out","r")
    ind=0

    for line in hand:
        line = line.rstrip()
        #[R00] Times: IO: 12.864861; Setup: 0.396950; Compute: 14.529508; MPI: 14.461607; Total: 27.980158;

        if re.search('^.*:.*:.*;$',line):
            number = re.findall('[-+]?([0-9]*\.[0-9]+|[0-9]+)',line)
            data[i,ind/120,ind%120,:]=np.array(number[1:])
            ind = ind+1


#http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html



batch_data=np.zeros((len(Size),4,num_file_to_merge))
batch_mean=np.zeros((4))
batch_std=np.zeros((4))
for i in range(len(Size)):
    #print(data[:,i,0:Pros[0],4].mean() , data[:,i,0:Pros[0],4].var() )
    #print(data[:,i,Pros[0]:Pros[1],4].mean(), data[:,i,Pros[0]:Pros[1],4].var() )
    #print(data[:,i,Pros[1]:Pros[2],4].mean(), data[:,i,Pros[1]:Pros[2],4].var() )
    #print(data[:,i,Pros[2]:Pros[3],4].mean(), data[:,i,Pros[2]:Pros[3],4].var() )
    for k in range(num_file_to_merge):
        batch_data[i,:,k]=np.array([data[k,i,0:Pros[0],plot_ind].mean(),
                                    data[k,i,Pros[0]:Pros[1],plot_ind].mean(),
                                    data[k,i,Pros[1]:Pros[2],plot_ind].mean(),
                                    data[k,i,Pros[2]:Pros[3],plot_ind].mean() ])
    for j in range(4):
        batch_mean[j] = batch_data[i,j,:].mean()
        batch_std[j] = batch_data[i,j,:].std()
    print batch_mean
    print batch_std
#    ax.set_title("system size: "+str(Size[i]))


