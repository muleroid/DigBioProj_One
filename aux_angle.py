import math
import numpy as np

#Finds angle ABC
def getAngle(ax,ay,az,bx,by,bz,cx,cy,cz, deg):
	print "ORDER IS IMPORTANT (getAngle). A,B,C"
	a = (bx - cx, by - cy)
	b = (bx - ax, by - ay)
	theta = np.arccos(np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b)))
	if(deg == True):
		return np.degrees(theta)
	return theta

#finds the distance between two points
def mag(mx,my,mz,nx,ny,nz):
	dsq=(mx-nx)*(mx-nx)+(my-ny)*(my-ny)+(mz-nz)*(mz-nz) 
	dist=math.sqrt(dsq)
	return math.sqrt(dist)


'''
def getAngle(ax,ay,az,bx,by,bz,cx,cy,cz, deg):
	#angle = arccos ( (av^2 + bv^2 - cv^2)/ -2avbv )
	av = mag(bx,by,bz,cx,cy,cz)
	bv = mag(ax,ay,az,cx,cy,cz)
	cv = mag(ax,ay,az,bx,by,bz)
	av2 = np.square(av)
	bv2 = np.square(bv)
	cv2 = np.square(cv)
	rad = np.arccos( (av2 + bv2 - cv2) / (-2*av*bv) )
	if(deg == True):
		return np.degrees(rad)
	return rad
'''