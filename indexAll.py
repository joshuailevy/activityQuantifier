#importing libraries/functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy import stats
from math import floor
from math import sqrt
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
import xlrd
import os
import sys
import csv


inputf = open('parameters.txt','r')
p1=[]
for line in inputf:
     name,value = line.split(':')
     p1.append(value.strip())
inputf.close()
dayType= int(p1[0])
initialZT=float(p1[1])
ZTon=float(p1[2])
ZToff=float(p1[3])
skew=int(p1[4])
dataType=int(p1[5])
hoursTotal=float(p1[6])
hoursBetween=float(p1[7])
def baseline_als(y, lam, p, niter=100):
  L = len(y)
  D = sparse.csc_matrix(np.diff(np.eye(L), 2))
  w = np.ones(L)
  for i in xrange(niter):
    W = sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = spsolve(Z, w*y)
    w = p * (y > z) + (1-p) * (y < z)
  return z

def individualExponential(t,coefs):
    return 0.5*(np.sign(t-coefs[2]) + 1)*coefs[0]*np.exp(-coefs[1]*(t-coefs[2])) 

def model(t,coefs):
    global dayType
    global noGauss1
    global noGauss2
    if noGauss1==0 and noGauss2==0:

    	if dayType==1:
       		return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2)) + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2)) \
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF)) \
       		+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2)) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2)) \
       		+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2)) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2)) 
    	else:
        	return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2)) + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2)) \
        	+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2)) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2)) \
        	+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2)) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2)) 
    
    elif noGauss1==1 and noGauss2==0:
    	if dayType==1:
       		return coefs[12] + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2)) \
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF)) \
       		+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2)) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2)) 
    	else:
        	return coefs[12] + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2)) \
        	+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2)) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2)) 

    elif noGauss1==0 and noGauss2==1:
    	if dayType==1:
       		return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2)) \
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF)) \
       		+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2)) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2))
    	else:
        	return coefs[12]+ coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2)) \
        	+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2)) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2))

    else:
    	if dayType==1:
       		return coefs[12] \
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF))
    	else:
        	return coefs[12]


def modelSkew(t,coefs):
    global dayType
    global noGauss1
    global noGauss2
    if noGauss1==0 and noGauss2==0:

    	if dayType==1:
       		return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-coefs[4])/(sqrt(2)*coefs[5])))\
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF)) \
       		+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t+24-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-24-coefs[1])/(sqrt(2)*coefs[2]))) \
       		+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t+24-coefs[4])/(sqrt(2)*coefs[5]))) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-24-coefs[4])/(sqrt(2)*coefs[5])))
    	else:
        	return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-coefs[4])/(sqrt(2)*coefs[5])))\
       		+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t+24-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-24-coefs[1])/(sqrt(2)*coefs[2]))) \
       		+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t+24-coefs[4])/(sqrt(2)*coefs[5]))) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-24-coefs[4])/(sqrt(2)*coefs[5])))
    	
    
    elif noGauss1==1 and noGauss2==0:
    	if dayType==1:
       		return coefs[12] + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-coefs[4])/(sqrt(2)*coefs[5])))\
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF)) \
       		+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t+24-coefs[4])/(sqrt(2)*coefs[5]))) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-24-coefs[4])/(sqrt(2)*coefs[5])))
       	else:
         	return coefs[12] + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-coefs[4])/(sqrt(2)*coefs[5])))\
       		+ coefs[3]*np.exp(-0.5*(((t+24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t+24-coefs[4])/(sqrt(2)*coefs[5]))) + coefs[3]*np.exp(-0.5*(((t-24-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-24-coefs[4])/(sqrt(2)*coefs[5])))
    	
    elif noGauss1==0 and noGauss2==1:
    	if dayType==1:
       		return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-coefs[1])/(sqrt(2)*coefs[2]))) \
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF)) \
       		+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t+24-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-24-coefs[1])/(sqrt(2)*coefs[2])))
    	else:
        	return coefs[12] + coefs[0]*np.exp(-0.5*(((t-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[3]*np.exp(-0.5*(((t-coefs[4])/coefs[5])**2))*(1+erf(coefs[14]*(t-coefs[4])/(sqrt(2)*coefs[5])))\
       		+ coefs[0]*np.exp(-0.5*(((t+24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t+24-coefs[1])/(sqrt(2)*coefs[2]))) + coefs[0]*np.exp(-0.5*(((t-24-coefs[1])/coefs[2])**2))*(1+erf(coefs[13]*(t-24-coefs[1])/(sqrt(2)*coefs[2])))

    else:
    	if dayType==1:
       		return coefs[12] \
       		+ 0.5*(np.sign(t-xON) + 1)*coefs[6]*np.exp(-coefs[7]*(t-xON)) + 0.5*(np.sign(t-xOFF) + 1)*coefs[9]*np.exp(-coefs[10]*(t-xOFF))
    	else:
        	return coefs[12]

    #+ coefs[6]*np.exp(-coefs[7]*(t-coefs[8])) + coefs[9]*np.exp(-coefs[10]*(t-coefs[11]))

#defining error function for our optimization problem 
def residuals(coefs, y, t):
    return (y - model(t,coefs))**2 +  100*(2 - np.sign(coefs[0]) -np.sign(coefs[3]))**2 + 10*(1-np.sign(coefs[12])) + 100*(coefs[12]>1.0) + 100*(coefs[6]<0)+ 100*(coefs[9]<0.5) + 100*(coefs[6]<0.5) + 100*(coefs[7]<0) + 100*(coefs[10]<0) +  0.01*(abs(coefs[2])>2)*abs(coefs[2])+ 100*(coefs[7]>1) + 100*(coefs[10]>1)#+ 0.1*(abs(coefs[2])>2)*abs(coefs[2])#+ 100*(1-np.sign(coefs[9]))**2 
    
#defining error function for our optimization problem 
def residualsSkew(coefs, y, t):
    return (y - modelSkew(t,coefs))**2 +  100*(2 - np.sign(coefs[0]) -np.sign(coefs[3]))**2  + 10*(1-np.sign(coefs[12]))  + 100*(coefs[12]>2.0) + 100*(coefs[6]<0.5)+ 100*(coefs[9]<0.5) +  100*(coefs[7]<0) + 100*(coefs[10]<0)+ 0.01*(abs(coefs[2])>2)*abs(coefs[2]) + 100*(coefs[7]>1) + 100*(coefs[10]>1)#+ 0#+ 0.1*(abs(coefs[2])>2)*abs(coefs[2])#+ 100*(coefs[12]>2.0) + 0.1*(abs(coefs[2])/abs(coefs[0]+0.1)>4)*abs(coefs[2])/abs(coefs[0]+0.1)# + 1*(abs(coefs[13]>1)+ abs(coefs[14])>1) 
    
def individualGaussian(t,iCoefs):
    return iCoefs[0]*np.exp(-0.5*(((t-iCoefs[1])/iCoefs[2])**2)) + iCoefs[0]*np.exp(-0.5*(((t+24-iCoefs[1])/iCoefs[2])**2)) + iCoefs[0]*np.exp(-0.5*(((t-24-iCoefs[1])/iCoefs[2])**2)) 

def individualSkewGaussian(t,iCoefs):
	return iCoefs[0]*np.exp(-0.5*(((t-iCoefs[1])/iCoefs[2])**2))*(1+erf(iCoefs[3]*(t-iCoefs[1])/(sqrt(2)*iCoefs[2]))) \
    + iCoefs[0]*np.exp(-0.5*(((t+24-iCoefs[1])/iCoefs[2])**2))*(1+erf(iCoefs[3]*(t+24-iCoefs[1])/(sqrt(2)*iCoefs[2]))) \
    + iCoefs[0]*np.exp(-0.5*(((t-24-iCoefs[1])/iCoefs[2])**2))*(1+erf(iCoefs[3]*(t-24-iCoefs[1])/(sqrt(2)*iCoefs[2])))

MI=[]
MIind = []
MIindSEM=[]
EI=[]
EIind = []
EIindSEM=[]
maxM= []
maxE=[]
peakLocM=[]
peakLocE=[]
pDist=[]
runName=[]
mOn=[]
mOff=[]
eOn=[]
eOff=[]
oldPhaseM=[]
oldPhaseE=[]

fileNum=0
noGauss1 = 0
noGauss2 = 0

time = np.linspace(0.0,hoursTotal,int(round(hoursTotal/hoursBetween)),endpoint=False)
Xs = time%24.0
x2 = np.linspace(0.0,24.0,int(round(24.0/hoursBetween)),endpoint=False)
xON = (ZTon)%24.0#(ZTon-initialZT)%24.0
xOFF = (ZToff)%24.0#(ZToff-initialZT)%24.0

if dataType == 1:
    for file in os.listdir("."):
        if file.endswith(".xls"):
            print file
            strippedfname = file.rsplit('.',1)[0]
            book = xlrd.open_workbook(file)
            #print book.nsheets
            #print book.sheet_names()
            sheet = book.sheet_by_index(0)
            fileNum = fileNum+1
            rowAvg = []
            MI1 = []
            EI1 = []
            avg6 =0.0
            avg3 = 0.0
            avg62 =0.0
            avg32 = 0.0
            in6avg=0
            in6avg2=0
            for i in xrange(56,sheet.nrows):
                avg=0.0;
                numAdded=0;
                #take second sheet. do averaging as before. 
                for j in xrange(0,sheet.ncols):
                    #print sheet.cell_type(rowx=i,colx=j)
                    if sheet.cell_type(rowx=i,colx=j):
                        avg = avg+sheet.cell_value(rowx=i,colx=j)
                        numAdded = numAdded+1
                rowAvg.append(avg/numAdded)

            Ys = np.array(rowAvg,dtype=float)
            for j in xrange(0,sheet.ncols): 
                for i in xrange(56,sheet.nrows): 
                    i1 = i-56 
                    tPoint = (xON -  float(i1)*hoursBetween)%(24.0)
                    tPoint2 = (xOFF -  float(i1)*hoursBetween)%(24.0)
                    #print tPoint, i1, (xON -  float(i1)*hoursBetween)
                    #compute morning anticipation for the individual flies
                    if ((tPoint <= 6) and tPoint > 0):
                        in6avg= in6avg +1
                        avg6 = avg6 + sheet.cell_value(rowx=i,colx=j)
                        #print avg6
                        if tPoint <= 3:
                            avg3 = avg3 + sheet.cell_value(rowx=i,colx=j)
                    if in6avg==(6/hoursBetween):
                        #avg3 = avg3#/(3.0/hoursBetween)
                        #avg6 = avg6/(6.0/hoursBetween)
                        if avg6 !=0.0:
                            bla = avg3/avg6
                            MI1.append(bla)
                        else:
                            MI1.append(0.0)
                        avg3=0.0
                        avg6=0.0
                        in6avg = 0
                 
                    #print tPoint, i1, (xON -  float(i1)*hoursBetween)f
                    #compute evening anticipation for the individual flies
                    if ((tPoint2 <= 6) and tPoint2 > 0):
                        in6avg2= in6avg2 +1
                        avg62 = avg62 + sheet.cell_value(rowx=i,colx=j)
                        if tPoint2 <= 3.0:
                            avg32 = avg32 + sheet.cell_value(rowx=i,colx=j)
                    if in6avg2==(6/hoursBetween):
                        #avg32 = avg32/(3.0/hoursBetween)
                        #avg62 = avg62/(6.0/hoursBetween)
                        if avg62 !=0.0:
                            blab = avg32/avg62
                            EI1.append(blab)
                        else:
                            EI1.append(0.0)
                        avg32=0.0
                        avg62=0.0
                        in6avg2 = 0
            
            print 'MI = ', np.mean(MI1), '+-', stats.sem(MI1)
            print 'EI = ', np.mean(EI1), '+-', stats.sem(EI1)
            #plt.hist(MI)
            #plt.show()
            #bline = baseline_als(Ys,1.0E4,0.01)
            #Ys = Ys-0.5
            #Ys = Ys-bline;#[0];
            #now average over days..
#plt.hist(EI1)
 #           plt.show()
            days = round(len(Xs)/48)
            w = np.zeros(48)
            for i in xrange(len(Xs)):
                ind = int(Xs[i]*2)
                w[ind] = w[ind] + Ys[i]

                #if Xs[i] == 0.0:
                #   days = days +
            w1 = w/days
            w = w/days


            ## calculate published phase measurement

            indON = np.argwhere(Xs==xON)
            indOFF = np.argwhere(Xs==xOFF)

           # print "on/off = ", indON, indOFF
            wEarly = w[0:indON[0][0]]
            wMiddle = w[indON[0][0]:indOFF[0][0]]
            overallVal = -1

            for i in xrange(len(wEarly)-3):
                val1 = wEarly[i+3]-wEarly[i]
                if val1 > overallVal:
                    ind11 = i
                    overallVal=val1

            oldMorningPhase = Xs[ind11]
            #print "time phase morning=", Xs[ind11]

            overallVal = -1
            ind11 = 0
            for i in xrange(len(wMiddle)-3):
                val1 = wMiddle[i+3]-wMiddle[i]
                if val1 > overallVal:
                    ind11 = i
                    overallVal=val1
           
            oldEveningPhase = Xs[ind11+indON[0][0]]
            #print "time phase evening=", Xs[ind11+indON[0][0]]


            #bline = baseline_als(w,1.0E5,0.01)
            #w=w-bline[0]
            print fileNum
            if fileNum==1:
                picking =1
            else:
                picking=0

            mVal =np.mean(w) - 0.5
            mR = 100000
           #x0 = [float(peakLocY[0]), float(peakLoc[0]), 1.0, float(peakLocY[1]),float(peakLoc[1]),1.0,2.0,2.0,xON,2.0,2.0,xOFF,0.3]#this can be modified... but should still give us conversion to the best fit
            #x0 = [float(peakLocY[0]), float(peakLoc[0]), 1.0, float(peakLocY[1]),float(peakLoc[1]),1.0,2.0,2.0,xON,2.0,2.0,xOFF,0.3]#this can be modified... but should still give us conversion to the best fit
            for i2 in range(-2,3,2):
                x0 = [1.0, xON-i2, 1.5, 1.0,xOFF,1.0,2.0,2.0,xON + .5,2.0,2.0,xOFF,mVal]#this can be modified... but should still give us conversion to the best fit
                #x0 = [1.0, float(peakLoc[0]), 1.0, 1.0,float(peakLoc[1]),1.0,5.0,4.0,float(peakLocSwitch[0]),5.0,4.0,float(peakLocSwitch[1])]#this can be modified... but should still give us conversion to the best fit
                res1 = least_squares(residuals,x0,loss='soft_l1',f_scale=1.0, args=(w,x2))#,verbose=1)
                if (i2==-2) or res.cost < mR:
                    res = res1

            x3 = np.linspace(x2[0],x2[len(x2)-1],10*len(x2))
            y1 = np.zeros(len(x3))
            y2 = np.zeros(len(x3))
            y3 = np.zeros(len(x3))
            y4 = np.zeros(len(x3))
            y5 = np.zeros(len(x3))
            #print res.x
            x= res.x

            if (x[2]<0.5 or x[2] > 5) or x[0] <0.5:
            	noGauss1 =1
            	res = least_squares(residuals,x0,loss='soft_l1',f_scale=1.0, args=(w,x2))
            	x= res.x
            	x[0] = 0
            	noGauss1=0

            if (x[5]<0.5 or x[5] > 5)  or x[3] <0.5:
            	noGauss2 =1
            	res = least_squares(residuals,x0,loss='soft_l1',f_scale=1.0, args=(w,x2))
            	x= res.x
            	x[3] = 0
            	noGauss2=0

            if ((x[5]<0.5 or x[5] > 5)  or x[3] <0.5) and ((x[2]<0.5 or x[2] >5) or x[0] <0.5):
            	noGauss2 =1
            	noGauss1 =1
            	res = least_squares(residuals,x0,loss='soft_l1',f_scale=1.0, args=(w,x2))
            	x= res.x
            	x[0] = 0
            	x[3] = 0
            	noGauss1=0
            	noGauss2=0

            if skew ==1:
            	x0Skew = np.append(x,(0.0,0.0))#[1.0, xON-1, 1.0, 1.0,xOFF-1,1.0,2.0,2.0,xON,2.0,2.0,xOFF,mVal,0.0,0.0]
            	if x[0]==0:
            		noGauss1=1
            		print "no morning peak"
            	if x[3]==0:
            		noGauss2=1
            		print "no evening peak"

            	resSkew = least_squares(residualsSkew,x0Skew,loss='soft_l1',f_scale=1, args=(w,x2))
            	x = resSkew.x
                print x
                if noGauss1 ==1:
                    x[0]=0
                if noGauss2 ==1:
                    x[3]=0

            	noGauss1=0
            	noGauss2=0

            for i in xrange(len(x3)):

	            if skew==0:
	                y1[i] =individualGaussian(x3[i],x[0:3])+x[12]
	                y2[i] =individualGaussian(x3[i],x[3:6])+x[12]
	                y3[i] =individualExponential(x3[i],x[6:9])+x[12]
	                y4[i] =individualExponential(x3[i],x[9:12])+x[12]
	                y5[i] = model(x3[i],x)
	            else:
	            	first = np.append(x[0:3],x[13])
	            	second = np.append(x[3:6],x[14])
	            	y1[i] =individualSkewGaussian(x3[i],first)+x[12]
	                y2[i] =individualSkewGaussian(x3[i],second)+x[12]
	                y3[i] =individualExponential(x3[i],x[6:9])+x[12]
	                y4[i] =individualExponential(x3[i],x[9:12])+x[12]
	                y5[i] = modelSkew(x3[i],x)
            #plt.plot(x2,y1)
            #plt.plot(x2,y2)
            #plt.plot(x2,y3)
            #plt.plot(x2,y4)
            
            #phase calculation by finding location of max and min of derivative
            y1diff = np.diff(y1)
            y1p = np.argmax(y1diff)
            y1n = np.argmin(y1diff)
            #print "new morning phase= ",x3[y1p], "offset = ", x3[y1n]

            y2diff = np.diff(y2)
            y2p = np.argmax(y2diff)
            y2n = np.argmin(y2diff)
            #print "new evening phase=", x3[y2p],"offset = ", x3[y2n]
            ax3 = plt.axes()
            plt.plot(2*x3,y5,color='r')
            #plt.plot(x2,y1)
            #plt.plot(x2,y2)
            #plt.bar(range(len(w)),w,align='center')
            tRange = range(len(w))
            #tRange = [i/2 for i in tRange1]
            plt.bar(tRange[0:2*int(round(xON))],w[0:2*int(round(xON))],align='center',color='k')
            plt.bar(tRange[2*int(round(xON)):2*int(round(xOFF))],w[2*int(round(xON)):2*int(round(xOFF))],align='center',color='w')
            plt.bar(tRange[2*int(round(xOFF)):len(w)],w[2*int(round(xOFF)):len(w)],align='center',color='k')
            #plt.xticks(np.arange(0,2*max(x3),4.0))#2*x3,x3,4.0)
            ax3.set_xticks(2*x2)
            #print [abs(i - round(i)) if abs(i - round(i))<0.0001 else abs(i - round(i)) for i in x3]
            ax3.set_xticklabels([('ZT'+str(int(i+ initialZT)%24) if i%4==0 else ' ') for i in x2])
            #plt.show()
            plt.savefig(strippedfname+'_model.png')
            plt.clf()


            ax = plt.axes()
            ax2 = plt.axes()
            if skew==0:
                p1x = x[1]
                p1y = x[0]+x[12]
                p2x = x[4]
                p2y = x[3]+x[12]
            else:
                ind1 = np.argmax(y1)
                ind2 = np.argmax(y2)
                p1x = x3[ind1]
                p1y = y1[ind1]
                p2x = x3[ind2]
                p2y = y2[ind2]

            if p1y-x[12]>0.01:
                #print p1y-x[12]
                ax.arrow(2*p1x,p1y+0.2,0,-0.1,head_width=0.3, head_length=0.1, fc='r', ec='r')
                ax.arrow(2*x3[y1p],y1[y1p]+0.2,0,-0.1,head_width=0.3, head_length=0.1, fc='b', ec='b')
                ax.arrow(2*x3[y1n],y1[y1n]+0.2,0,-0.1,head_width=0.3, head_length=0.1, fc='g', ec='g')

            if p2y-x[12]>0.01:
                #print p2y-x[12]
                ax2.arrow(2*p2x,p2y+0.2,0,-0.1,head_width=0.3, head_length=0.1, fc='r', ec='r')
                ax2.arrow(2*x3[y2p],y2[y2p]+0.2,0,-0.1,head_width=0.3, head_length=0.1, fc='b', ec='b')
                ax2.arrow(2*x3[y2n],y2[y2n]+0.2,0,-0.1,head_width=0.3, head_length=0.1, fc='g', ec='g')
            #plt.plot(x2,w,"o")
            
            ax.set_xticks(2*x2)
            #print [abs(i - round(i)) if abs(i - round(i))<0.0001 else abs(i - round(i)) for i in x3]
            ax.set_xticklabels([('ZT'+str(int(i + initialZT)%24) if i%4==0 else ' ') for i in x2])
            #print xON
            plt.bar(tRange[0:2*int(round(xON))],w[0:2*int(round(xON))],align='center',color='k')
            plt.bar(tRange[2*int(round(xON)):2*int(round(xOFF))],w[2*int(round(xON)):2*int(round(xOFF))],align='center',color='w')
            plt.bar(tRange[2*int(round(xOFF)):len(w)],w[2*int(round(xOFF)):len(w)],align='center',color='k')
            plt.plot(2*x3,y1,color='r')
            plt.plot(2*x3,y2,color='r')
            plt.savefig(strippedfname+'.png')
            plt.clf()
            #sys.exit()
            #output information for further inspection
    
            if x[0]!=0 and x[3]!=0:
                peakLocM.append(p1x)
                peakLocE.append(p2x)
                maxM.append(p1y)
                maxE.append(p2y)
                pDist.append(p2x-p1x)
                mOn.append(x3[y1p])
                mOff.append(x3[y1n])
                eOn.append(x3[y2p])
                eOff.append(x3[y2n])
                oldPhaseM.append(oldMorningPhase)
                oldPhaseE.append(oldEveningPhase)
            elif x[0]==0 and x[3]!=0:
                peakLocM.append(-1)
                peakLocE.append(p2x)
                maxM.append(0)
                maxE.append(p2y)
                pDist.append(-1)
                mOn.append(-1)
                mOff.append(-1)
                eOn.append(x3[y2p])
                eOff.append(x3[y2n])
                oldPhaseM.append(oldMorningPhase)
                oldPhaseE.append(oldEveningPhase)
            elif x[3]==0 and x[0]!=0:
                peakLocM.append(p1x)
                peakLocE.append(-1)
                maxM.append(p1y)
                maxE.append(0)
                pDist.append(-1)
                mOn.append(x3[y1p])
                mOff.append(x3[y1n])
                eOn.append(-1)
                eOff.append(-1)
                oldPhaseM.append(oldMorningPhase)
                oldPhaseE.append(oldEveningPhase)
            else:
                peakLocM.append(-1)
                peakLocE.append(-1)
                maxM.append(0)
                maxE.append(0)
                pDist.append(-1)
                mOn.append(-1)
                mOff.append(-1)
                eOn.append(-1)
                eOff.append(-1)
                oldPhaseM.append(oldMorningPhase)
                oldPhaseE.append(oldEveningPhase)

            #calculate morning and evening indices
            runName.append(strippedfname)
            mSw= np.argwhere(abs(x2-float(xON))<0.01)[0][0]
            eSw= np.argwhere(abs(x2-float(xOFF))<0.01)[0][0]
            last6m = np.sum(w1[mSw-12:mSw])
            last3m = np.sum(w1[mSw-6:mSw])
            last6e = np.sum(w1[eSw-12:eSw])
            last3e = np.sum(w1[eSw-6:eSw])
            MI.append(last3m/last6m)
            MIind.append(np.mean(MI1))
            MIindSEM.append(stats.sem(MI1))
            EIind.append(np.mean(EI1))
            EIindSEM.append(stats.sem(EI1))
            EI.append(last3e/last6e)

            book.release_resources()
            del book
            	#print book.nsheets
            	#print book.sheet_names()
                #print file#(os.path.join(subdir, file))

    fout = open('dataCalculated.csv','wb')
    try:
        writer=csv.writer(fout)
        writer.writerow(('File name','Morning Index','Ind. MI','Ind MI SEM','Evening Index','Ind. EI','Ind EI SEM','Peak loc. (morning)','Morning ant. height','Peak loc. (evening)','Evening ant. height','Dist. between peaks','morning onset','morning offset','evening onset','evening offset','published phase (morn.)','published phase (evening)'))
        for i in range(0,fileNum):
            writer.writerow((runName[i],MI[i],MIind[i],MIindSEM[i],EI[i],EIind[i],EIindSEM[i],peakLocM[i],maxM[i],peakLocE[i],maxE[i],pDist[i],mOn[i],mOff[i],eOn[i],eOff[i],oldPhaseM[i],oldPhaseE[i]))
    finally:
        fout.close()
else:
    print "This data input type is not yet supported. \n"
    
'''
fname = 'winter5min.csv'
strippedfname = fname.rsplit('.',1)[0]
f = open(fname,'r')
#extracting the data
xVals = []
yVals = []
i = 0#9#21
#reading the file
for line in f:
	line = line.strip()
	columns = line.split()
	yVals.append(columns[0])
	xVals.append(i%24)
	i=i+(1.0/12.0)
#plt.plot(xVals,yVals)
#Xs = np.array(xVals,dtype=float)
Ys = np.array(yVals,dtype=float)
Xs = np.array(xVals,dtype=float)
#plt.plot(Xs,Ys,"x")
'''
#print days
#plt.plot(xrange(len(w)),w,'o')
#plt.plot(Xs,Ys,'o')
#plt.show()

#x2 = np.array(Xs[xrange(288)],dtype=float)

'''


'''
