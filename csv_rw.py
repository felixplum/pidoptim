import csv
import numpy as np
import time
def read_model_csv(path):
    
    with open(path, 'r', newline='') as csvfile:
        file = csv.reader(csvfile, delimiter=';')
        data =[r for r in file]
        

        data =data[2:-1]

        temp=np.zeros((len(data),6)).astype(np.float)
        timee=list()
        for i in range(0,len(data)):
            temp2=data[i]
            temp3=time.strptime(temp2[0], "%d.%m.%Y %H:%M")
     
            timee.append(temp3[3]*60*60+temp3[4]*60)            
            
            temp[i,0]=float(temp2[2].replace(',', '.'))
            temp[i,1]=float(temp2[3].replace(',', '.' ))
            temp[i,2]=float(temp2[4].replace(',', '.'))
            temp[i,3]=float(temp2[5].replace(',', '.'))
            temp[i,4]=float(temp2[6].replace(',', '.' ))
            temp[i,5]=float(temp2[7].replace(',', '.'))
        
        #data=np.array(data).astype(np.float)
        
        t0=timee[0]
        timee=[item -t0 for item in timee]
        timee = np.array(timee,np.float)
      
        for i in range(0,int(timee[len(timee)-1]+1),60):
            tmp = timee==i
            tmp2 = timee[tmp]
            tmp3 = np.zeros(len(timee))
            timee[tmp]=np.linspace(i,i+60,num=len(tmp2),endpoint=False)#np.linspace(i,i+1,num=len(tmp2),endpoint=False)
            
        t_x=timee#np.linspace(0,timee[len(timee)-1],num=timee[len(timee)-1]/0.1)    
        
        t_x= t_x[t_x>10000]
        #t_x= t_x[t_x<19000]

        r = np.interp(t_x, timee, temp[:,0]) 
        y = np.interp(t_x, timee, temp[:,1]) 
        u = np.interp(t_x, timee, temp[:,2]) 
     

        csvfile.close()
        return t_x,u,y,r
    
    


