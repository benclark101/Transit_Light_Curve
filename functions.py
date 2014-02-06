from userinput import *

#AlignImages - A function that calculates the misalignment of a set of images in pixels up to +/- 50px
#Only works in low populated fields (becomes inaccurate if many stars in image)
#Simply requires location of FITS file folder
#Outputs an array to file of correction values in (x,y) pixels
def AlignImages(fitsLoc):
    #Get list of sources    
    Sources = []
    SourceList = csv.reader(open('sources.txt'), delimiter=',')
    j = 0
    for row in SourceList:
        X, Y = row
        x = float(X)
        y = float(Y)
        Sources.append([x,y])
        j = j + 1
    fileNames = os.listdir(fitsLoc)
    filenumber = 1
    for name in fileNames:
        print name
        filenumber = filenumber + 1
        filename = str(fitsLoc) + str(name)
        image  = pyfits.getdata(filename,0)
        TotalCountOld = 0
        for correctionX in range(-50,50):
            for correctionY in range(-50,50):
                TotalCount = 0
                for i in range(0,j):
                    InNew = 0
                    Y = Sources[i][0]
                    X = Sources[i][1]
                    X = X + correctionX
                    Y = Y + correctionY
                    for sq in range(-5,5):
                        InNew = InNew + float(image[X+sq,Y+sq])
                    TotalCount = TotalCount + InNew
                if TotalCount > TotalCountOld:
                    BestCorX = correctionX
                    BestCorY = correctionY
                    TotalCountOld = TotalCount
        fout = open("AlignArray.txt", "a")
        fout.write(str(name)+str(":")+str(BestCorX)+str(":")+str(BestCorY))
        fout.write("\n")            
        print BestCorX, BestCorY

#AperturePhotometry - Function That performs aperture photometry with a circular aperture
#Inputs: x,y - coordinates of source (pixels)
#r - radius of desired aperture (px)
#image - FITS imnage imported into python by  pyfits.getdata("FITS FILE LOCATION",0)
#Returns the count (total pixel value) within the region
def AperturePhotometry(x, y, r, image):
    #Create array same size as image
    R_arr = numpy.ones((image.shape[0],image.shape[1]))
    R_arr = R_arr * 4000
    #Calculate max and min points away from x and y coords (just above radius)
    refPixMaxX = int(x + r + .5)
    refPixMinX = int(x - r - .5)
    refPixMaxY = int(y + r + .5)
    refPixMinY = int(y - r - .5) 
    #Calculate distance away from center x,y point
    for i in range(refPixMinY, refPixMaxY):
        for j in range (refPixMinX, refPixMaxX):
            R_arr[i,j] = ((i-y)**2+(j-x)**2)**.5
    #select only pixels within given radius and sum        
    selectAper = numpy.where( R_arr <=  rMax  )
    sourceCount = image[selectAper].sum()
    Intensity = sourceCount
    return Intensity
           
#Center - Ensures that apertures are centered over brightest pixel (ie. apertures centered over stars)
#Inputs: x,y - coordinates of source (px)
#r - radius of desired aperture (px)
#image - FITS imnage imported into python by  pyfits.getdata("FITS FILE LOCATION",0)
#Returns centred coordinates where  aperture will be placed
def Center(x,y,r,image):
    refPixMaxX = int(x + r + 1)
    refPixMinX = int(x - r - 1)
    refPixMaxY = int(y + r + 1)
    refPixMinY = int(y - r - 1) 
    InOld = 0
    for i in range(refPixMinY, refPixMaxY):
        for j in range (refPixMinX, refPixMaxX):
            InNew = int(image[i,j])
            if InNew > InOld:
                Xcen = j
                Ycen = i
                InOld = InNew
    return (Xcen,Ycen)
    
#GetIntensity - Retrieve background subtracted Intensities of all sources located in Sources.txt
#Inputs - filename:location of file
#a - Iteration factor from "Main" code
#corX, corY - correction factors produced by alignment function (imported for the text file output by the code)
#Outputs the values of BG subtracted intensity to file
def GetIntensity(fileName,a, corX, corY,bgX,bgY):
    #Define Images
    image  = pyfits.getdata(fileName,0)
    bg = AperturePhotometry(bgX+corX,bgY+corY , rMax, image)
    #Get time of observation from header    
    FitsFile = pyfits.open(fileName)
    prihdr = FitsFile[0].header
    Date = prihdr['DATE-OBS']
    dte = Date.split("T")
    time = dte[1].split(":")
    hours = time[0]
    mins = time[1]
    secs = time[2]
    MPM = (int(hours) * 60) + int(mins) + (float(secs)/60)
    #Write observation time to file
    fout = open(IntenOutNm, "a")
    fout.write(str(MPM))
    fout.write(":")
    #Open list of sources
    SourceList = csv.reader(open('Sources.txt'), delimiter=',')
    for row in SourceList:
        X, Y = row
        x = float(X)
        y = float(Y)
        x = x + corX
        y = y + corY
        #Recenter source coords and calculate intensities
        x, y = Center(x, y, rMax, image)
        inten = AperturePhotometry(x,y, rMax, image)
        Intensity = inten - bg
        #Write values to file
        fout.write(str(x))
        fout.write(":")
        fout.write(str(y))
        fout.write(":")
        fout.write(str(Intensity))
        fout.write(":")
    if a == 0:
        fout = open(IntenOutNm, "a")
        fout.write("\n")
        fout.close
        
#<<<<<<< HEAD
#ErApPhotom - Performs aperture photometry for the error calculations (same as AperturePhotometry() above except the code calculates the standard deviation in the aperture too)
#=======
#ErApPhotom - Aperture photometry for error calculations
#same as AperturePhotometry() function above except this function also returns the standard deviation in the aperture
#>>>>>>> Minor Changes
def ErApPhotom(x, y, r, image):
    R_arr = numpy.ones((image.shape[0],image.shape[1]))
    R_arr = R_arr * 4000
    refPixMaxX = int(x + r + .5)
    refPixMinX = int(x - r - .5)
    refPixMaxY = int(y + r + .5)
    refPixMinY = int(y - r - .5)
    OutArr = []
    for i in range(refPixMinY, refPixMaxY):
        for j in range (refPixMinX, refPixMaxX):
            R_arr[i,j] = ((i-y)**2+(j-x)**2)**.5
            OutArr.append(image[i,j])
    SdCount = numpy.std(OutArr)
    selectAper = numpy.where( R_arr <=  rMax  )
    sourceCount = image[selectAper].sum()
    Intensity = sourceCount
    return Intensity, SdCount
    
#Errors - Same as GetIntensity() function but with error calculations
#Outputs error values to file
def Errors(fileName,a,corX, corY,bgX,bgY):
    image  = pyfits.getdata(fileName,0)
    bg,erBg = ErApPhotom(bgX+corX,bgY+corY , rMax, image)
    FitsFile = pyfits.open(fileName)
    prihdr = FitsFile[0].header
    Date = prihdr['DATE-OBS']
    dte = Date.split("T")
    time = dte[1].split(":")
    hours = time[0]
    mins = time[1]
    secs = time[2]
    MPM = (int(hours) * 60) + int(mins) + (float(secs)/60)
    fout = open("error.txt", "a")  
    fout.write(str(MPM))
    fout.write(":")
    SourceList = csv.reader(open('sources.txt'), delimiter=',')
    for row in SourceList:
        X, Y = row
        x = float(X)
        y = float(Y)
        x = x + corX
        y = y + corY
        x, y = Center(x, y, rMax, image)
        inten, StEr = ErApPhotom(x,y, rMax, image)
        Intensity = numpy.sqrt(bg + inten + numpy.square(StEr))
        fout.write(str(x))
        fout.write(":")
        fout.write(str(y))
        fout.write(":")
        fout.write(str(Intensity))
        fout.write(":")
    if a == 0:
        fout = open("error.txt", "a")
        fout.write("\n")
        fout.close        
        
#paramCalc - Calculates exoplanet parameters
#Requires inputs from the Radius() function below as well as ingress and egress times
#+++++ Need more accurate way of retrieving ingress and egress times +++++
#Outputs an array containing: Impact Parameter, Period, Semi-Major Axis, Inclination Angle, Mass and Density
def paramCalc(iS,iE,eS,eE,ratio,ErRatio, Rplanet, ErRplanet):
    #Constants
    G = 6.67E-11
    ExoParams = []
    #Calculate tf, tt and Delta
    ingressStart = float(iS) * 60
    ingressEnd = float(iE) *60
    egressStart = float(eS) * 60
    egressEnd = float(eE) * 60
    ratio = Rplanet / Rstar
    D = ratio**2
    ErD = D * numpy.sqrt((ErRatio/ratio)**2 +(ErRatio/ratio)**2 )
    rootD = numpy.sqrt(D)
    ErRtD = rootD * numpy.sqrt(((0.5*ErD)/D)**2)
    Tf = float(egressStart - ingressEnd)
    ErTf = float(200) #(estimated error in seconds)
    Tt = float(egressEnd - ingressStart)
    ErTt = float(200)
    Tratio = Tf/Tt
    ErTratio = float(Tratio * numpy.sqrt((ErTf/Tf)**2 + (ErTt/Tt)**2 ))
    #Impact Parameter calculation
    b = numpy.sqrt( (((1-rootD)**2 - Tratio**2) * (1+rootD)**2 ) / (1 - Tratio**2))
    ErB = b * numpy.sqrt(  2*((ErRtD/rootD)**2) + 2*((ErTratio/Tratio)**2 )    )
    ExoParams.append(b)
    #Period of orbit
    Period = (Mstar/(Rstar**3)) * ((G*3.14159265358979)/32) * ( ((Tt**2 - Tf**2)**1.5) / D**0.75 )
    PeriodD = Period / (3600*24) #days
    TermAEr = (ErMstar/Mstar**2) + ((3*ErRstar)/Rstar)**2
    TermBEr = (3*ErTf/Tf)**2 + (3*ErTt/Tt)**2 + (0.75*ErD/D)**2
    ErPeriod = Period * numpy.sqrt(TermAEr + TermBEr)
    ErPeriodD = ErPeriod /(3600*24)
    print "Period: ", PeriodD, " +/- ", ErPeriodD, "Days"
    ExoParams.append(Period)
    #Semi-major Axis
    a = sp.cbrt(((Period**2) * G * Mstar) / (4*((3.14159265358979)**2)))
    aAu = a / 1.5E11 #Au
    ErA = a * numpy.sqrt((0.666*float(ErPeriod)/float(Period))**2 + (0.333*float(ErMstar)/float(Mstar))**2)
    ErAau = ErA / 1.5E11
    print "Semi-major Axis: ", aAu, " +/- ", ErAau, "Au"
    ExoParams.append(a)
    #Inclination Angle
    incl = numpy.arccos((b*Rstar)/a)
    inclD = numpy.rad2deg(incl)
    ErIncl = inclD * numpy.sqrt((ErB/b)**2 + (ErRstar/Rstar)**2 + (ErA/a)**2)
    print "Inclination Angle: ", inclD, " +/- ", ErIncl, "Degrees"
    ExoParams.append(incl)
    #Mass of planet
    #Mrv = 1.404 * 1.898E27
    #ErMrv = 0.099 * 1.898E27
    Mass = Mrv / numpy.sin(incl)
    ErMass = Mass * numpy.sqrt((ErMrv/Mrv)**2 + (numpy.deg2rad(ErIncl)*numpy.cos(incl))**2 )
    MassJ = Mass / 1.898E27
    ErMassJ = ErMass / 1.898E27
    print "Mass: ", MassJ," +/- ", ErMassJ," Mj"
    ExoParams.append(Mass)
    #Density of planet
    Density = Mass / (1.33*3.14159265358979*(Rplanet**3))
    ErDensity = Density * numpy.sqrt( (ErMass/Mass)**2 + (0.33*ErRplanet/Rplanet)**2)
    Densitygcm = Density/1000
    ErDengcm = ErDensity/1000
    print "Density: ", Densitygcm, " +/- ", ErDengcm, "gcm^-3"
    ExoParams.append(Density)
    print "************************************************"
    #Exoparams: b, P, a, incl, mass, density
    return ExoParams

#Radius - Calculates the radius of the exoplanet using the dip in luminoisity and radius of host star
def Radius(mtop, mbot, stdtop, stdbot,Ttimes):
    DeltaF = mtop-mbot
    ErrDF = numpy.sqrt(numpy.square(stdtop*stdtop) + numpy.square(stdbot*stdbot))
    Foff = mtop
    ErrFO = stdbot
    LumRatio = (DeltaF/Foff)
    ELumRatio = LumRatio * numpy.sqrt((ErrDF/DeltaF) + (ErrFO/Foff))
    #Rstar = 1.57 * 6.955E8
    Rplanet = numpy.sqrt(LumRatio * (Rstar * Rstar))
    #ERstar = 0.07 * 6.955E8
    ErrRP = Rplanet * numpy.sqrt(numpy.square(ELumRatio/LumRatio) + numpy.square(ERstar / Rstar))
    RplanJ = Rplanet / 71492000
    ERPj = ErrRP / 71492000
    RplanJ = round(RplanJ, 3)
    ERPj = round(ERPj, 3)
    print 
    print "************************************************"
    print "Radius is: ", '{0:.3f}'.format(RplanJ), "Jupiter Radii +/- ", '{0:.3f}'.format(ERPj)
    fout = open("params.txt", "w")
    fout.write(str("Calculated Radius: ") + str(RplanJ) + str(" Jupiter Radii +/- ") + str(ERPj) + str("\n"))
    fout.write(str("LumRatio: ") + str(LumRatio) + str("\n"))
    TiA = Ttimes[0]
    TiB = Ttimes[1]
    TeA = Ttimes[2]
    TeB = Ttimes[3]
    ParamArray = paramCalc(TiA,TiB,TeA,TeB,LumRatio,ELumRatio,Rplanet,ErrRP)
    fout.write(str(ParamArray))
    return ParamArray

#find - Find matching time values between the data and model arrays    
def find(j, intensity, i, FluxArr):  
    IntNew = []
    FluxNew = []
    accuracy = 0.1 #Accuracy to with which arrays are matched (in seconds)
    F = 0
    for b in j:
        a = 0
        for c in i:
            if b >= c - accuracy and b<= c + accuracy:
                IntNew.append(intensity[F])
                FluxNew.append(FluxArr[a])
                break
            a = a + 1
        F = F + 1
    #print len(IntNew)
    if len(IntNew) < (0.75*len(j)):
        print "WARNING - Matches too low - Increase Z points in modelPlot()* or change accuracy in find()"
    IntNP = numpy.asarray(IntNew)
    FlNP = numpy.asarray(FluxNew)
    value = numpy.sum((IntNP - FlNP)**2)
    return value

         
#inten - Intensity equation I(r) from Mandel & Agol 2002
def inten(r,cn):
    IrSum = 0
    a = float(cn[0] * (1-((1-(r*r))**0.25 )))
    b = float(cn[1] * (1-((1-(r*r))**0.5 )))
    c = float(cn[2] * (1-((1-(r*r))**0.75 )))
    d = float(cn[3] * (1-((1-(r*r))**1 )))
    IrSum = a + b + c + d
    Ir = float((1. - IrSum) * 2. * r)
    return Ir
    
#omega - Omega equation from Mandel & Agol 2002
def omega(cn):
    Omg = 0
    coeffs = []
    c0 = 1 - cn[0] - cn[1] - cn[2] - cn[3]
    coeffs.append(c0)
    for i in range(0,4):
        coeffs.append(cn[i])
    n = 0
    for CN in coeffs:
        a = float(CN) / float(n+4)
        n = n + 1
        Omg = Omg + a
    return Omg

#z2t - Convert from 'z' to time
def z2t(tmid, Fullz):
    TimeArray = []
    P = PeriodVal
    a = SemiMjrAxis
    i = InclAngle
    for z in Fullz:
        if z > 0:
            t = tmid + ( ((P*Rstar) / (2*3.1415926535*a)) * (z**2 - (numpy.cos(i))**2)**0.5)
        if z < 0:
            t = tmid + ( ((P*Rstar) / (2*3.1415926535*a)) * -(z**2 - (numpy.cos(i))**2)**0.5)
        TimeArray.append(t)
    return TimeArray



#modelPlot - Function that creates and fits model    
def modelPlot(zarr, intensity,IntenErrP, params, parameters):
    p = params[0]
    Tmid = params[1]
    iter8 = 0
    #retrieve exoplanet parmeters
    global PeriodVal
    global SemiMjrAxis
    global InclAngle
    SemiMjrAxis = parameters[2]
    InclAngle = parameters[3]
    PeriodVal = parameters[1]
    #Host star limb darkening coefficients (from spectral info):
    CNarray = [[0.4857,0.2088,0.2177,-0.1753]]
    Coefficients = str("0.4857,0.2088,0.2177,-0.1753")
    for cn in CNarray:
        iter8 = iter8 + 1
        IntenZ = []
        Zarray = []
        FluxArr = []
        #*Change POINTS below - more 'z' points gives better accuracy when performing modelling but takes more time
        POINTS = 6000
        for z in numpy.linspace(0.00000001,2.55,POINTS):
            if z < 1 - p:
                LowerLim = z - p
                UpperLim = z + p
                a = (z - p) ** 2
                IntensityIntgr,null = intg.quad(lambda x: inten(x,cn), LowerLim, UpperLim)
                IZ = (IntensityIntgr) / (4*z*p)
                IntenZ.append(IZ)
                Zarray.append(z)
                Om = omega(cn)
                Flux = 1 - (((1/(4*Om)) * IZ) * p ** 2)
                FluxArr.append(Flux)
            elif z < 1 + p and z > 1 - p:
                LowerLim = z - p
                UpperLim = 1
                a = (z - p) ** 2
                IntensityIntgr,null = intg.quad(lambda x: inten(x,cn), LowerLim, UpperLim)
                IZ =  IntensityIntgr / (1-a)
                IntenZ.append(IZ)
                Zarray.append(z)
                Om = omega(cn)
                Flux = 1 - ((IZ / (4*Om*numpy.pi)) * ((p*p*numpy.arccos((z-1)/p)) - (z-1)*numpy.sqrt((p*p) - ((z-1)**2))))
                FluxArr.append(Flux)
            else:
                FluxArr.append(1)
                Zarray.append(z)
        Zarray2 =  []
        Zarray3 =  []
        for z in Zarray:
            Zarray2.append(-z)
        for z in Zarray:
            Zarray3.append(z)
        TimeArray = []
        OldTimeArr = []
        Fullz =  Zarray2 + Zarray3
        OldLstSq = 100000000
        OldT = 0
        FluxArr = FluxArr + FluxArr
        Tmin = Tmid - 1000
        Tmax = Tmid + 1000
        step = 100
        for looper in [0,1,2]:
            if looper == 0:
                print "Stage 1/3"
            elif looper == 1:
                print "Stage 2/3"
            elif looper == 2:
                print "Stage 3/3"
            TimeArr = z2t(Tmin, Fullz)
            To = numpy.asarray(TimeArr)
            FullLSarray = []
            for Tmiddle in range(Tmin,Tmax,step):
                TimeArray = To + (Tmiddle - Tmin)
                LeastSqVal = find(zarr, intensity, TimeArray, FluxArr)
                if LeastSqVal < 0.01:
                    LeastSqVal = LeastSqVal + 0.04
                FullLSarray.append((Tmiddle, LeastSqVal))
            valu0 = 10E9
            Tbest = 0
            for valu in FullLSarray:
                if valu[1] < valu0:
                    valu0 = valu[1]
                    Tbest = valu[0]
            if looper == 0:
                Tmin = Tbest - 100
                Tmax = Tbest + 100
                step = 10
                looper = 1
            elif looper == 1:
                Tmin = Tbest - 10
                Tmax = Tbest + 10
                step = 1
                looper = 2
            elif looper == 2:
                break
        TimeArray = z2t(Tbest, Fullz)
        plt.plot(TimeArray, FluxArr, '-', linewidth=2, label= 'Model with Coeffs: %s' % Coefficients)
        plt.plot(zarr, intensity, 'ro')
        plt.errorbar(zarr, intensity, yerr=IntenErrP ,linestyle="None", marker="None", color='grey')
    return (Tbest, valu0)
