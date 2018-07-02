import csv
import numpy as np

def read_model_csv(path):
    
    with open(path, 'r', newline='') as csvfile:
        file = csv.reader(csvfile, delimiter=',')
        data =[r for r in file]
        
        data=np.array(data[1:-1]).astype(np.float)
        t = data[:,0]
        u = data[:,1]
        y = data[:,2] 
        r = data[:,3]         
        csvfile.close()
        return t,u,y,r
    
    


