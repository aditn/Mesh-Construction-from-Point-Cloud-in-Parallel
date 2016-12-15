from random import random
import math

def dist(p1,p2):
    (x1,y1,z1)=p1
    (x2,y2,z2)=p2
    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**0.5

def makeSphereMesh(c,r,filename,numPoints=10000):
    """c: tuple of 3, r: int
    randomly samples points on sphere's surface"""
    points=[]
    for i in xrange(numPoints):
        theta = math.pi*random()
        phi = (2*math.pi*random())-math.pi
        x = r*math.sin(theta)*math.cos(phi)
        y = r*math.sin(theta)*math.sin(phi)
        z = r*math.cos(theta)
        points += [(c[0]+x,c[1]+y,c[2]+z)]

    f = open(filename,"w")
    for p in points:
        f.write("v %.3f %.3f %.3f\n" % p)
    f.close()

def isSphere(c,r,points):
    for p in points:
        if(abs(dist(c,p)-r)>0.001):
            return False
    return True

def getSquarePointsRandom(pos,negs,numPoints):
    w=pos[0]-negs[0]
    h=pos[1]-negs[1]
    d=pos[2]-negs[2]
    points = []
    for i in xrange(numPoints):
        points += [(negs[0]+random()*w,negs[1]+random()*h,negs[2]+random()*d)]

def getSquarePoints(pos,negs,numPoints):
    uniformps = int(math.floor(math.sqrt(numPoints)))-1
    w=pos[0]-negs[0]
    h=pos[1]-negs[1]
    d=pos[2]-negs[2]
    points = []
    if(w==0): #yz
        for i in xrange(uniformps+1):
            for j in xrange(uniformps+1):
                points += [(negs[0],negs[1]+float(i)/uniformps*h,negs[2]+float(j)/uniformps*d)]
    elif(h==0): #xz
        for i in xrange(uniformps+1):
            for j in xrange(uniformps+1):
                points += [(negs[0]+float(i)/uniformps*w,negs[1],negs[2]+float(j)/uniformps*d)]
    else: #xy
        for i in xrange(uniformps+1):
            for j in xrange(uniformps+1):
                points += [(negs[0]+float(i)/uniformps*w,negs[1]+float(j)/uniformps*h,negs[2])]
    return points

def makeCubeMesh(c,s,filename,numPoints=10000):
    """uniformally samples points on cube's surface"""
    points = []
    pps = numPoints/6 #points per side
    blb = (c[0]-s*0.5,c[1]-s*0.5,c[2]-s*0.5) #bottomleftback
    blf = (c[0]-s*0.5,c[1]-s*0.5,c[2]+s*0.5) #bottomleftback
    brb = (c[0]+s*0.5,c[1]-s*0.5,c[2]-s*0.5) #bottomleftback
    brf = (c[0]+s*0.5,c[1]-s*0.5,c[2]+s*0.5) #bottomleftback
    tlb = (c[0]-s*0.5,c[1]+s*0.5,c[2]-s*0.5) #bottomleftback
    tlf = (c[0]-s*0.5,c[1]+s*0.5,c[2]+s*0.5) #bottomleftback
    trb = (c[0]+s*0.5,c[1]+s*0.5,c[2]-s*0.5) #bottomleftback
    trf = (c[0]+s*0.5,c[1]+s*0.5,c[2]+s*0.5) #bottomleftback
    
    #xy planes
    points += getSquarePoints(blb,trb,pps)
    points += getSquarePoints(blf,trf,pps)
    #yz planes
    points += getSquarePoints(blb,tlf,pps)
    points += getSquarePoints(brb,trf,pps)
    #xz planes
    points += getSquarePoints(blb,brf,pps)
    points += getSquarePoints(tlb,trf,pps)

    f = open(filename,"w")
    for p in points:
        f.write("v %.3f %.3f %.3f\n" % p)
    f.close()
