from userinput import *
import functions as fn

#Create (/overwrite) new text file
fout = open(IntenOutNm, "w")
erout = open(ErrOutNm, "w")   
print "Processing Files:"
#Get list of filenames
fileNames = os.listdir(fitsLoc)
#Show image of first FITS file to allow for source locating
print "Open sources.txt and populate with X-Y coordinates of 10 star locations with the first being the target star - note background location"
filename = str(fitsLoc) + str(fileNames[0])
image  = pyfits.getdata(filename,0)
plt.figure(1)
plt.imshow(image)
#if stars not visible in image, change value of 1000 below (try multiples of 100 above or below until sources are visible)
plt.clim(0,1000)
plt.show()
global bgX
global bgY
bgX = input("Enter X background coordinate: ")
bgY = input("Enter Y background coordinate: ")
#Align all images
fn.AlignImages(fitsLoc)
#***********Get intensities for all sources in each image****************
filenumber = 1
for i in fileNames:
    print i
    filename = str(fitsLoc) + str(i)
    ArrayAlign = csv.reader(open('AlignArray.txt'), delimiter=':')
    for row in ArrayAlign:
        FileNameA, CorX, CorY = row
        if FileNameA == i:
            corX = float(CorY)
            corY = float(CorX)
            print i
    fn.GetIntensity(filename,a, corX, corY,bgX,bgY)
    fn.Errors(filename,a, corX, corY,bgX,bgY)
    filenumber = filenumber + 1
#Plot Light Curve and retrieve radius
FileNumb = filenumber
IntVar = 'Int'
ErVar = 'Err'
sourcenumb = 10
#Create arrays based on number of sources
charSET = []
alphabet = ['S','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','T','U','V','W','X','Y','Z']
for iterate in range(0,sourcenumb):
    charSET.append(alphabet[iterate])
for char in charSET:
    vars()[str(IntVar) + str(char)] = numpy.zeros(FileNumb)
    vars()[str(ErVar) + str(char)] = numpy.zeros(FileNumb)
data = csv.reader(open('data.txt'), delimiter=':')
dataErr1 = csv.reader(open('error.txt'), delimiter=':')
Time = numpy.zeros(FileNumb)
i = 0
#Read Intensity Data into 48arrays
for row in data:
    MPM, xS, yS, IS, xA, yA, IA, xB, yB, IB, xC, yC, IC, xD, yD, ID, xE, yE, IE, xF, yF, IF, xG, yG, IG, xH, yH, IH, xI, yI, Ii, null = row
    IntS[i] = float(IS)
    IntA[i] = float(IA)
    IntB[i] = float(IB)
    IntC[i] = float(IC)
    IntD[i] = float(ID)
    IntE[i] = float(IE)
    IntF[i] = float(IF)
    IntG[i] = float(IG)
    IntH[i] = float(IH)
    Time[i] = float(MPM)
    i = i + 1 
#Find Observation Start Time    
LowestTime = 1000000000
for t in range(0,FileNumb):
    if Time[t] < LowestTime:
        LowestTime = Time[t]
Time = Time - LowestTime
i = 0
#Read Error values into arrays
for row in dataErr1:
    MPM, xS, yS, IS, xA, yA, IA, xB, yB, IB, xC, yC, IC, xD, yD, ID, xE, yE, IE, xF, yF, IF, xG, yG, IG, xH, yH, IH, xI, yI, Ii, null = row
    ErrS[i] = float(IS)
    ErrA[i] = float(IA)
    ErrB[i] = float(IB)
    ErrC[i] = float(IC)
    ErrD[i] = float(ID)
    ErrE[i] = float(IE)
    ErrF[i] = float(IF)
    ErrG[i] = float(IG)
    ErrH[i] = float(IH)   
    i = i + 1
#Find best calibrator - lowest standard deviation
IA = numpy.mean([numpy.std(IntA/IntB),numpy.std(IntA/IntC),numpy.std(IntA/IntD),numpy.std(IntA/IntE)])
IB = numpy.mean([numpy.std(IntB/IntC),numpy.std(IntB/IntD),numpy.std(IntB/IntE),numpy.std(IntB/IntA)])
IC = numpy.mean([numpy.std(IntC/IntD),numpy.std(IntC/IntE),numpy.std(IntC/IntA),numpy.std(IntC/IntB)])
ID = numpy.mean([numpy.std(IntD/IntE),numpy.std(IntD/IntA),numpy.std(IntD/IntB),numpy.std(IntD/IntC)])
IE = numpy.mean([numpy.std(IntE/IntA),numpy.std(IntE/IntB),numpy.std(IntE/IntC),numpy.std(IntE/IntD)])
#Find lowest mean - get best calibrators
MeanArray = [(IA, "A"), (IB, "B"), (IC,  "C"), (ID, "D"), (IE, "E")]
MeanArray.sort()
#Create Super Calibrator
Numb1, IntLowest1 = MeanArray[0]
Numb2, IntLowest2 = MeanArray[1]
Numb3, IntLowest3 = MeanArray[2]
superCal = (vars()[str('Int') + str(IntLowest1)] + vars()[str('Int') + str(IntLowest2)] + vars()[str('Int') + str(IntLowest3)]) / 3
superCalEr = numpy.sqrt((numpy.square(vars()[str('Err') + str(IntLowest1)]) + numpy.square(vars()[str('Err') + str(IntLowest2)]) + numpy.square(vars()[str('Err') + str(IntLowest3)])) / 3)
#Compute intensity ratios with error
Intensity = IntS / superCal
IntenErrP = Intensity * numpy.sqrt(numpy.square(ErrS /IntS) + numpy.square(superCalEr/superCal) )
I = str('ISuperCal')
#Calculate mean and standard deviation of light curve
print "Enter Ingress and egress times"
plt.plot(Time, Intensity, 'ro')
plt.show()
TingressA = input("Ingress start: ")
TingressB = input("Ingress end: ")
TegressA = input("Egress start: ")
TegressB = input("Egress end: ")
TransitTiming = [TingressA,TingressB,TegressA,TegressB]
ingressA, ingressB, egressA, egressB = 0,0,0,0
for T in Time:
    if T <= TingressA:
        ingressA = ingressA + 1
    if T <= TingressB:
        ingressB = ingressB + 1
    if T <= TegressA:
        egressA = egressA + 1
    if T <= TegressB:
        egressB = egressB + 1
Itop = []
Ibot = []
for i in range (0,ingressA):
    Itop.append(Intensity[i])
for i in range (egressB,FileNumb):
    Itop.append(Intensity[i])
for i in range (ingressB,egressA):
    Ibot.append(Intensity[i])
mtop = numpy.mean(Itop)
stdTop = numpy.std(Itop)
mbot = numpy.mean(Ibot)
stdBot = numpy.std(Ibot)
#******************Calculate Exoplanet Parameters********************
parameters = fn.Radius(mtop,mbot,stdTop,stdBot,TransitTiming)    
#Normalise Intensity Ratio    
Intensity = Intensity / mtop 
IntenErrP = IntenErrP / mtop
#Divide throughout by mean to get 1 mean value
stdTop = stdTop / mtop
mbot = mbot / mtop
stdBot = stdBot / mtop
mtop = 1
#********************Plot Light Curve********************
plt.figure(1, figsize=(10, 6))
plt.plot(Time, Intensity, 'ro')
plt.errorbar(Time, Intensity, yerr=IntenErrP ,linestyle="None", marker="None", color='grey')
plt.axis(xmin = -20,xmax = 360)
plt.axhline(y=mtop, xmin=0, xmax=0.21, hold=None,linewidth=2, color='k')
plt.axhline(y=mtop+stdTop, xmin=0, xmax=0.21, hold=None, color='r')
plt.axhline(y=mtop-stdTop, xmin=0, xmax=0.21, hold=None, color='r')
plt.axhline(y=mtop, xmin=0.60, xmax=1, hold=None,linewidth=2, color='k')
plt.axhline(y=mtop+stdTop, xmin=0.60, xmax=1, hold=None, color='r')
plt.axhline(y=mtop-stdTop, xmin=0.60, xmax=1, hold=None, color='r')
plt.axhline(y=mbot, xmin=0.20, xmax=0.65, hold=None,linewidth=2, color='k')
plt.axhline(y=mbot+stdBot, xmin=0.20, xmax=0.65, hold=None, color='b')
plt.axhline(y=mbot-stdBot, xmin=0.20, xmax=0.65, hold=None, color='b')
#Format and save plot to file
plt.axis(xmin=-10,ymin = 0.975,ymax = 1.01)
plt.title('Light Curve for Wasp-12b')
plt.xlabel('Time /"mins after observation start"')
plt.ylabel(str('Normalised Intensity Ratio (IS/') + str(I)+ str(') / Mean'))
plt.savefig("Graph.png", format = 'png')
plt.show()
plt.clf()
#Write arrays to file
fout = open("Zvals.txt", "w")
Zarr = Time * 60
for Zs in Zarr:
    fout.write(str(Zs) + str(',0')+ str('\n'))
fout = open("Intensities.txt", "w")
for Zs in Intensity:
    fout.write(str(Zs) + str(',0') + str('\n'))
fout = open("IntenErrors.txt", "w")
lupe = 0
for Zs in IntenErrP:
    lupe = lupe + 1
    fout.write(str(Zs) + str(',0'))
    if lupe == len(IntenErrP):
        DoNothing = True
    else:
        fout.write(str('\n'))
fout.close()
#********************Create and Fit Model to data********************
print "Creating and Fitting Model"
i = 0
New = csv.reader(open('Zvals.txt'), delimiter=',')
for row in New:
    i = i + 1
j = 0
New2 = csv.reader(open('Intensities.txt'), delimiter=',')
for row in New2:
    j = j + 1  
Zarray = numpy.ones(i)
IntArray = numpy.ones(j)
i = 0
SourceList = csv.reader(open('Zvals.txt'), delimiter=',')
for row in SourceList:
    Zr, r = row
    Zr = float(Zr)
    Zarray[i] = Zr
    i = i + 1
i = 0    
SourceList2 = csv.reader(open('Intensities.txt'), delimiter=',')
for row in SourceList2:
    Int, r = row
    Int = float(Int)
    IntArray[i] = Int
    i = i + 1
Znew = Zarray
BestArray = []
for Pnumb in numpy.linspace(0.119,0.121,3):
    Tguess = 8300
    params = [Pnumb,Tguess]    
    BestT, BestVal = fn.modelPlot(Zarray,IntArray,0, params, parameters)
    BestArray.append((Pnumb,BestT,BestVal))
print BestArray
iold = 10000000
for i in BestArray:
    if i[2] < iold:
        iold = i[2]
        BestTMid = i[1]
        BestP = i[0]
print BestTMid
print BestP
print iold
plt.show()   
#**************Perform Monte Carlo Simulations********************
print "Performing Monte Carlo Simulations (This may take a while...)"
#++++++++Efficiency of algorithm needs improving++++++++
for looper44 in range(0,200): 
    print "MAIN : ", looper44
    i = 0
    New = csv.reader(open('Zvals.txt'), delimiter=',')
    for row in New:
        i = i + 1
    j = 0
    New2 = csv.reader(open('Intensities.txt'), delimiter=',')
    for row in New2:
        j = j + 1    
    Zarray = numpy.ones(i)
    IntArray = numpy.ones(j)
    Err = numpy.ones(j)
    i = 0
    SourceList = csv.reader(open('Zvals.txt'), delimiter=',')
    for row in SourceList:
        Zr, r = row
        Zr = float(Zr)
        Zarray[i] = Zr
        i = i + 1 
    i = 0    
    SourceList2 = csv.reader(open('Intensities.txt'), delimiter=',')
    for row in SourceList2:
        Int, r = row
        Int = float(Int)
        IntArray[i] = Int
        i = i + 1
    i = 0    
    SourceList2 = csv.reader(open('IntenErrors.txt'), delimiter=',')
    for row in SourceList2:
        Errr, r = row
        Errr = float(Errr)
        Err[i] = Errr
        i = i + 1
    Inew = []
    iter8 = 0
    for Xi in IntArray:
        Error = Err[iter8]
        Error = int(Error * 10000)
        Inew.append(Xi + (float(random.randint(-Error,Error)) / 10000))
        iter8 = iter8 + 1
    BestArray = []
    for Pnumb in numpy.linspace(0.11,0.13,21):
        Tguess = 8300
        params = [Pnumb,Tguess]       
        BestT, BestVal = fn.modelPlot(Zarray,Inew,0, params, parameters)
        BestArray.append((Pnumb,BestT,BestVal))
    print BestArray
    iold = 10000000
    for i in BestArray:
        if i[2] < iold:
            iold = i[2]
            BestTMid = i[1]
            BestP = i[0]
    fout = open("BestValAr.txt", "a")
    fout.write(str(BestTMid)+ str(",")+ str(BestP)+  str(",")+str(iold)+str(",") +str("\n"))
    print "Saved"
#Plot frequency and find mean and std
Inew = []
New = csv.reader(open('BestValArr.txt'), delimiter=',')
i = 0
for row in New:
    a,b,c,d=row
    Inew.append(float(a))
    i = i+ 1
print len(Inew)    
plt.hist(Inew,bins=20)
print numpy.mean(Inew)
print numpy.std(Inew)
plt.title('Monte Carlo Simulation - Transit Midpoint Distribution')
plt.xlabel('Transit Midpoint')
plt.ylabel(str('Number'))
plt.axvline(x=numpy.mean(Inew), color='c')
plt.axvline(x=numpy.mean(Inew)+numpy.std(Inew), color='r')
plt.axvline(x=numpy.mean(Inew)-numpy.std(Inew), color='r')
plt.show()  
#********************Create O-C Plot********************
#+++++++++NEEDS TO BE AUTOMATED - VALUES CURRENTLY INPUT MANUALLY+++++++++++++
input("O-C plot currently requires manual editing - please change values in main.py before continuing")
ArchiveData = csv.reader(open('C:/Users/ben/Desktop/wasp12ocdata.csv'), delimiter=',')
epocha = []
OCa = []
for row in ArchiveData:
    epoch, OC = row
    epocha.append(float(epoch))
    OCa.append(float(OC))
plt.plot(epocha, OCa, 'ro')
plt.axhline(y=0)
plt.axis(ymin = -0.015,ymax = 0.02)
plt.axis(ymin = -0.015,ymax = 0.011)
#plt.axis(xmin=1200, xmax=1400)
plt.plot(1301,-0.007777, 'bs')
plt.errorbar(1301, -0.0077, yerr=0.003067 ,linestyle="None", marker="None", color='grey')
plt.plot(1302,-0.00903, 'bs')
plt.errorbar(1302,-0.00903, yerr=0.002396 ,linestyle="None", marker="None", color='grey')
plt.plot(1327,-0.007777, 'bs')
plt.errorbar(1327,-0.007777, yerr=0.002697 ,linestyle="None", marker="None", color='grey')
plt.title('O-C plot for Wasp-12b')
plt.xlabel('Epoch')
plt.ylabel(str('Observed-Calculated /days'))
plt.show()   
#****************Attempt to fit basic sin curve to plot****************
sinar = []
epoch2 = []
OCb = []
epochb = []
i = 0
for epochs in epocha:
    if epochs > 1200 and epochs <1400:
        epochb.append(epochs)
        OCb.append(OCa[i])
        i = i + 1      
OCb = numpy.asarray(OCb)
epochb = numpy.asarray(epochb)
fitfunc = lambda p, x: p[0]*numpy.sin(p[1]*x) + p[2]
errfunc = lambda p, x, y: fitfunc(p, x) - y 
p0 = [0.005, 0.1,1] 
p1, success = optimize.leastsq(errfunc, p0, args=(epochb, OCb),maxfev=100000000)
print p1
for i in numpy.linspace(1200,1400,1000):
    sinv = p1[0]*numpy.sin(i*p1[1])+p1[2]
    sinar.append(sinv)
    epoch2.append(i)
plt.plot(epoch2,sinar)
sinar = []
epoch2 = []
OCa = numpy.asarray(OCa)
epocha = numpy.asarray(epocha)
fitfunc = lambda p, x: p[0]*numpy.sin(p[1]*x) + p[2]
errfunc = lambda p, x, y: fitfunc(p, x) - y 
p0 = [0.1, 0.1,1] 
p1, success = optimize.leastsq(errfunc, p0, args=(epocha, OCa))
print p1
for i in numpy.linspace(1200,1400,1000):
    sinv = p1[0]*numpy.sin(i*p1[1])+p1[2]
    sinar.append(sinv)
    epoch2.append(i)

plt.plot(epoch2,sinar)

plt.show() 
