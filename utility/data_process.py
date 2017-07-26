__author__="Xiang Liu"
from numpy import random
import numpy as np
import pandas as pd
from simpleFunctions import *
from numpy import linalg
from matplotlib import pyplot as plt

##here how to process the real data(traits)
a = pd.read_csv("../toyData/alzheimer/traits.csv",header=None,keep_default_na=0)  #keep_default_na=0
train_data = a.values
##
print train_data.shape

temp=train_data[:,0]
temp[temp=='NA']=0.
train_data[:,0]=temp
###########0------------##############
temp=train_data[:,1]
temp[temp=='NA']=0.
train_data[:,1]=temp
###########1------------##############
temp=train_data[:,2]
temp[temp=='NA']=0.
train_data[:,2]=temp
###########2------------###############
temp=train_data[:,3]
temp[temp=='NA']=0.
train_data[:,3]=temp
###########3------------###############
temp=train_data[:,4]
temp[temp=='NA']=0.
train_data[:,4]=temp
###########4------------###############
temp=train_data[:,5]
temp[temp=='NA']=0.
train_data[:,5]=temp
###########5------------###############
temp=train_data[:,6]
temp[temp=='NA']=0.
train_data[:,6]=temp
###########6------------###############
temp=train_data[:,7]
temp[temp=='NA']=0.
train_data[:,7]=temp
###########7------------###############
temp=train_data[:,8]
temp[temp=='NA']=0.
train_data[:,8]=temp
###########8------------###############
temp=train_data[:,9]
temp[temp=='NA']=0.
train_data[:,9]=temp
###########9------------###############
temp=train_data[:,10]
temp[temp=='NA']=0.
train_data[:,10]=temp
###########10------------###############
temp=train_data[:,11]
temp[temp=='NA']=0.
train_data[:,11]=temp
###########11------------###############
temp=train_data[:,12]
temp[temp=='NA']=0.
train_data[:,12]=temp
###########12------------###############
temp=train_data[:,13]
temp[temp=='NA']=0.
train_data[:,13]=temp
###########13------------###############
temp=train_data[:,14]
temp[temp=='NA']=0.
train_data[:,14]=temp
###########14------------###############
temp=train_data[:,15]
temp[temp=='NA']=0.
train_data[:,15]=temp
###########15------------###############
temp=train_data[:,16]
temp[temp=='NA']=0.
train_data[:,16]=temp
###########16------------###############
temp=train_data[:,17]
temp[temp=='NA']=0.
train_data[:,17]=temp
###########17------------###############
temp=train_data[:,18]
temp[temp=='NA']=0.
train_data[:,18]=temp
###########18------------###############
temp=train_data[:,19]
temp[temp=='NA']=0.
train_data[:,19]=temp
###########19------------###############
temp=train_data[:,20]
temp[temp=='NA']=0.
train_data[:,20]=temp
###########20------------###############
temp=train_data[:,21]
temp[temp=='NA']=0.
train_data[:,21]=temp
###########21------------###############
temp=train_data[:,22]
temp[temp=='NA']=0.
train_data[:,22]=temp
###########22------------###############
temp=train_data[:,23]
temp[temp=='NA']=0.
train_data[:,23]=temp
###########23------------###############
temp=train_data[:,24]
temp[temp=='NA']=0.
train_data[:,24]=temp
###########24------------###############
temp=train_data[:,25]
temp[temp=='NA']=0.
train_data[:,25]=temp
###########25------------###############
temp=train_data[:,26]
temp[temp=='NA']=0.
train_data[:,26]=temp
###########26------------###############
temp=train_data[:,27]
temp[temp=='NA']='no'
temp[temp=='no']=0.
temp[temp=='yes']=1.
train_data[:,27]=temp
###########27------------###############
temp=train_data[:,28]
temp[temp=='NA']=1.    #1.*56/108
temp[temp=='F']=0.    #52
temp[temp=='M']=1.    #56
train_data[:,28]=temp
###########28------------###############
temp=train_data[:,29]
temp[temp=='LNV']=0.
temp[temp=='Dry-Ice']=1.
train_data[:,29]=temp
###########29------------###############
########30
temp=train_data[:,31]
temp[temp=='control']=0.
temp[temp=='affected']=1.
train_data[:,31]=temp
###########31------------###############
temp=train_data[:,32]
temp[temp=='control']=0.
temp[temp=='affected']=1.
train_data[:,32]=temp
###########32------------###############
temp=train_data[:,33]
temp[temp=='Unknown']=0.
temp[temp=='Vonsattel']=2. #63
temp[temp=='VonSattel']=2. #63
temp[temp=='Lathi']=1.     #105
temp[temp=='Hedreen']=0.   #361
train_data[:,33]=temp
###########33------------###############
temp=train_data[:,34]
temp[temp=='HD 4']=0.2
temp[temp=='HD 3']=0.5
temp[temp=='HD 2']=0.8
temp[temp=='Control']=1.
temp[temp=='AD']=2.
temp[temp=='AD/Braak 6']=3.1
temp[temp=='AD/Braak 5']=3.2
temp[temp=='AD/Braak 4']=3.3
temp[temp=='AD/Braak 3']=3.4
temp[temp=='AD/Braak 2']=3.5
temp[temp=='AD/Braak 5/6']=3.6
temp[temp=='AH1']=4.1
temp[temp=='AH2']=4.2
temp[temp=='AH3']=4.3
temp[temp=='AH']=4.4
train_data[:,34]=temp
###########34------------###############
temp=train_data[:,35]
temp[temp=='NA']=0. #(5.*2+1.*135+2.*12+81+40)/(197+135+61)
temp[temp=='E2/E2']=5.  #2
temp[temp=='E3/E3']=0.  #197
temp[temp=='E3/E4']=1.  #135
temp[temp=='E2/E4']=2.  #12
temp[temp=='E2/E3']=3.  #27
temp[temp=='E4/E4']=4.  #20
train_data[:,35]=temp
###########35------------###############
for i in range(0,540):
    train_data[i,36]=float(train_data[i,36])
train_data[:,36]=normalize(train_data[:,36])
#########36
################################the last lines are greater than 1
temp=train_data[:,37]
temp[temp=='NA']=17.4021666667  #mean
for i in range(0,540):
   temp[i]=float(temp[i])
train_data[:,37]=temp
train_data[:,37]=normalize(train_data[:,37])
 #########37########
temp=train_data[:,38]
temp[temp=='NA']=46.1111111111
temp[temp=='NO DATA']=46.1111111111
for i in range(0,540):
    temp[i]=float(temp[i])
temp=normalize(temp)
train_data[:,38]=temp
########38###########
temp=train_data[:,39]
temp[temp=='NA']=18.9090909091
temp[temp=='NO DATA']=18.9090909091
for i in range(0,540):
    temp[i]=float(temp[i])
temp=normalize(temp)
train_data[:,39]=temp
##########39
temp=train_data[:,40]
temp[temp=='NA']=40.8275862069
for i in range(0,540):
    temp[i]=float(temp[i])
temp=normalize(temp)
train_data[:,40]=temp
##########40

temp=train_data[:,41]
temp[temp=='NA']=6.37480679406
for i in range(0,540):
    temp[i]=float(temp[i])
temp=normalize(temp)
train_data[:,41]=temp
##########41
train_data[train_data=='5.3b']=5.3
for i in range(0,540):
    train_data[i,42]=float(train_data[i,42])
train_data[:,42]=normalize(train_data[:,42])
for i in range(0,train_data.shape[0]):
    for j in range(0,train_data.shape[1]):
        train_data[i][j]=float(train_data[i][j])

# for i in range(0,train_data.shape[1]):
#     print "-------------- ",i," ------------------ "
#     print train_data[:,i]
np.savetxt('../toyData/after_process_2.csv', train_data, delimiter=',')
train_data=normalize(train_data)

np.savetxt('../toyData/after_process.csv', train_data, delimiter=',')




#here  how to process the real data(snps)
# fp=open("../toyData/alzheimer/snps.csv", "r");
# X=np.array([])
# time=0
# for eachline in fp:
#     x = []
#     ls = 0
#     time=time+1
#     print "times: ",time
#     for i in eachline:
#         if ls<=555090:
#             if i is  ',':
#                 pass
#             elif  i is ' ,':
#                 pass
#             elif i is ' ':
#                 pass
#             else:
#                 ls=ls+1
#                 i_1=float(i)
#                 x.append(i_1)
#         else:
#             pass
#     x=np.array(x)
#     x.reshape((1,555091))
#     if time==1 :
#         X=x
#     else:
#         X=np.row_stack((X, x))
#     #print x
# # X = np.array(X)
#
# # print "ok~ "
# X.reshape((540,555091))
# # X = np.load("../toyData/alzheimer/snps.npy")
# np.save("../toyData/alzheimer/snps.npy", X)

