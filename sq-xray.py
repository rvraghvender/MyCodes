import numpy as np

print(" Sequence O - Te - Tl ") 

c1 = float(input("Enter the value of c1 : "))  # Stoichiometric fractinal concentration
c2 = float(input("Enter the value of c2 : "))
c3 = float(input("Enter the value of c3 : "))

ci = [c1, c2, c3]

################# Form Factor Calculation #################                     
                                                                                
                                                                                
def form_factor(q):                                                             
    global Form_Factor                                                          
    Form_Factor = 0.0                                                           
    for ik in range(4):                                                         
        Form_Factor = Form_Factor  + ai[ik] * np.exp( -bi[ik] * ( q / (4*np.pi) )**2 )
    # print(Form_Factor + ci)                                                   
    Form_Factor = Form_Factor + c                                               
    return Form_Factor                                                          
                                                                                
                                                                                
def Te_form_factor(q):                                                          
    global ai                                                                   
    global bi                                                                   
    global c                                                                    
#    ai = np.array([6.660302, 6.940756, 19.847015, 1.557175,17.802427])                           
#    bi = np.array([33.031654, 0.025750, 5.065547, 84.101616,0.487660])                         
#    c  = -0.806668                                                                  
    ai = np.array([19.9644,19.0138, 6.14487, 2.5239])                           
    bi = np.array([4.81742,0.420885, 28.5284, 70.8403])                         
    c  = 4.352                                                                  
    form_factor(q)                                                              
    return Form_Factor                                                          
                                                                                
                                                                                
def O_form_factor(q):                                                           
    global ai                                                                   
    global bi                                                                   
    global c                                                                    
#    ai = np.array([ 2.960427, 2.508818, 0.637853, 0.722838,1.142756])                             
#    bi = np.array([14.182259, 5.936858, 0.112726, 34.958481,0.390240])                           
#    c  = 0.027014                                                                 
    ai = np.array([ 3.0485, 2.2868, 1.5463, 0.867])                             
    bi = np.array([13.2771, 5.7011, 0.3239, 32.9089])                           
    c  = 0.2508                                                                 
    form_factor(q)                                                              
    return Form_Factor                                                          
                                                                                 
def Tl_form_factor(q):                                                           
    global ai                                                                   
    global bi                                                                   
    global c                                                                    
#    ai = np.array([ 16.630795,19.386616,32.808571,1.747191,6.356862])     # Atomic factor constants for Tl1+                        
#    bi = np.array([0.110704,7.181401,1.119730,90.660263, 26.014978])                           
#    c  = 4.066939                                                                 
    ai = np.array([ 21.3985, 20.4723, 18.7478, 6.82847])     # Atomic factor constants for Tl1+                        
    bi = np.array([1.4711,0.517394, 7.43463, 28.8482])                           
    c  = 12.5258                                                                 
#    ai = np.array([ 27.5446, 19.1584, 15.538, 5.52593])     # Atomic factor constant for Tl                        
#    bi = np.array([0.65515, 8.70751, 1.96347, 45.8149])                           
#    c  = 13.1746                                                                 
    form_factor(q)                                                              
    return Form_Factor                                                          
                                                                                 
def Ti_form_factor(q):                                                           
    global ai                                                                   
    global bi                                                                   
    global c                                                                    
    ai = np.array([ 9.7595, 7.3558, 1.6991, 1.9021])     # Atomic factor constants for Tl1+                        
    bi = np.array([ 7.8508, 0.5   , 35.6338, 116.105])                           
    c  = 1.2807                                                                 
    form_factor(q)                                                              
    return Form_Factor                                                          

#############################################################  
#############################################################  

f1 = open("fz_O_O.dat",'r')
f2 = open("fz_O_Te.dat",'r')
f3 = open("fz_O_Ti.dat",'r')
f4 = open("fz_Te_Te.dat",'r')
f5 = open("fz_Te_Ti.dat",'r')
f6 = open("fz_Ti_Ti.dat",'r')

fw = open("sq-xr.dat",'w')

fq = [0,0,0]

while True:

    line1 = f1.readline().split()
    line2 = f2.readline().split()
    line3 = f3.readline().split()
    line4 = f4.readline().split()
    line5 = f5.readline().split()
    line6 = f6.readline().split()

    if len(line1) == 0 : break
    if line1[0] == '#' : continue

    sum_denominator = 0.0
    sum_numerator = 0.0

    fz_ij = [ float(line1[1]), float(line2[1]), float(line3[1]) , float(line4[1]) , float(line5[1]) , float(line6[1]) ]

    fq[0] = O_form_factor(float(line1[0]))                                     
    fq[1] = Te_form_factor(float(line1[0]))  
    fq[2] = Tl_form_factor(float(line1[0]))  

    k = 0

    for i in range(3):
        sum_denominator = sum_denominator + ci[i] * fq[i]**2

        for j in range(i,3):
            sum_numerator = sum_numerator + ci[i] * ci[j] * fq[i] * fq[j] * ( fz_ij[k] - 1 )
            k += 1

    sq = ( sum_numerator / sum_denominator ) + 1
    fw.write("	{}  	 {} \n".format(line1[0],sq))

fw.close()
