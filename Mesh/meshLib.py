import numpy as np

tol = 1.0E-10
Nround = 10

def setTol(tolValue):
    global tol 
    global Nround
    tol = tolValue
    Nround = -int(np.log10(tol))
    pass

class Vertex():
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = np.round(x, Nround)
        self.y = np.round(y, Nround)
        self.z = np.round(z, Nround)
        pass

def destVertex(v1, v2):
    return np.sqrt((v1.x-v2.x)**2+(v1.y-v2.y)**2+(v1.z-v2.z)**2)

def isEqualVertex(v1, v2):
    global tol
    if (np.abs(v1.x-v2.x)<tol) and\
        (np.abs(v1.y-v2.y)<tol) and\
        (np.abs(v1.z-v2.z)<tol):
        return True
    else:
        return False

class Line():
    def __init__(self, v1, v2, dL, portIndex):
        self.v1 = v1
        self.v2 = v2
        if dL>self.length():
            self.dL = self.length()
            pass
        else:
            self.dL = dL
            pass
        self.portIndex = portIndex
        pass
    def length(self):
        return destVertex(self.v1, self.v2)

class Segment():
    def __init__(self, v1, v2, portIndex, index):
        self.v1 = v1
        self.v2 = v2
        self.portIndex = portIndex
        self.index = index
        self.contiguous = []
        pass
    def __repr__(self):
        string = "v1 ({:21.14E}, {:21.14E}, {:21.14E})\n".format(self.v1.x, self.v1.y, self.v1.z)
        string+= "v2 ({:21.14E}, {:21.14E}, {:21.14E})\n".format(self.v2.x, self.v2.y, self.v2.z)
        return string

class Basis():
    def __init__(self, rm, rn, rp, portIndex):
        self.rm = rm
        self.rn = rn
        self.rp = rp
        self.portIndex = portIndex
        pass

class Shape():
    def __init__(self):
        self.basisList = []
        self.segmentList = []
        self.segmentIndex = 0
        pass
    def addLines(self, linesList):
        for line in linesList:
            n = np.ceil(line.length()/line.dL)
            try:
                assert(n>0)
            except:
                print("Error: Incorrect Dimensions!")
                exit(1)
            if np.mod(n, 2)!=0:
                n+=1
                pass
            line.dL = line.length()/n
            alpha = line.dL/line.length()
            for i in range(int(n)):
                v1 = Vertex(line.v1.x+(i+0)*alpha*(line.v2.x-line.v1.x),\
                            line.v1.y+(i+0)*alpha*(line.v2.y-line.v1.y),\
                            line.v1.z+(i+0)*alpha*(line.v2.z-line.v1.z))
                v2 = Vertex(line.v1.x+(i+1)*alpha*(line.v2.x-line.v1.x),\
                            line.v1.y+(i+1)*alpha*(line.v2.y-line.v1.y),\
                            line.v1.z+(i+1)*alpha*(line.v2.z-line.v1.z))
                if (line.portIndex>0 and (i+1)==int(n)/2) or\
                    line.portIndex>0 and (i+1)==1+int(n)/2:
                    newSegment = Segment(v1, v2, line.portIndex, self.segmentIndex)
                    pass
                else:
                    newSegment = Segment(v1, v2, 0, self.segmentIndex)
                    pass
                self.segmentIndex+=1
                self.segmentList.append(newSegment)
                continue
            continue
        pass
    def addSegments(self, segmentsList):
        for segment in segmentsList:
            segment.segmentIndex = self.segmentIndex
            self.segmentIndex+=1
            self.segmentList.append(segment)
            continue
        pass
    def getBasis(self):
        for segmentSource in self.segmentList:
            for segmentDest in self.segmentList:
                if segmentDest.index not in segmentSource.contiguous and \
                    segmentSource.index not in segmentDest.contiguous and \
                        segmentSource.index != segmentDest.index:
                    count = 0
                    if isEqualVertex(segmentSource.v1, segmentDest.v1):
                        if segmentSource.portIndex==segmentDest.portIndex:
                            newBasis = Basis(segmentSource.v2,\
                                             segmentSource.v1,\
                                             segmentDest.v2,\
                                             segmentSource.portIndex)
                            count+=1
                            pass
                        else:
                            newBasis = Basis(segmentSource.v2,\
                                             segmentSource.v1,\
                                             segmentDest.v2,\
                                             0)
                            count+=1
                            pass
                        segmentSource.contiguous.append(segmentDest.index)
                        segmentDest.contiguous.append(segmentSource.index)
                        try:
                            assert(count==1)
                        except:
                            print("Error: Incorrect Shape!")
                            exit(1)
                        self.basisList.append(newBasis)
                        pass
                    elif isEqualVertex(segmentSource.v1, segmentDest.v2):
                        if segmentSource.portIndex==segmentDest.portIndex:
                            newBasis = Basis(segmentSource.v2,\
                                             segmentSource.v1,\
                                             segmentDest.v1,\
                                             segmentSource.portIndex)
                            count+=1
                            pass
                        else:
                            newBasis = Basis(segmentSource.v2,\
                                             segmentSource.v1,\
                                             segmentDest.v1,\
                                             0)
                            count+=1
                            pass
                        segmentSource.contiguous.append(segmentDest.index)
                        segmentDest.contiguous.append(segmentSource.index)
                        try:
                            assert(count==1)
                        except:
                            print("Error: Incorrect Shape!")
                            exit(1)
                        self.basisList.append(newBasis)
                        pass
                    elif isEqualVertex(segmentSource.v2, segmentDest.v1):
                        if segmentSource.portIndex==segmentDest.portIndex:
                            newBasis = Basis(segmentSource.v1,\
                                             segmentSource.v2,\
                                             segmentDest.v2,\
                                             segmentSource.portIndex)
                            count+=1
                            pass
                        else:
                            newBasis = Basis(segmentSource.v1,\
                                             segmentSource.v2,\
                                             segmentDest.v2,\
                                             0)
                            count+=1
                            pass
                        segmentSource.contiguous.append(segmentDest.index)
                        segmentDest.contiguous.append(segmentSource.index)
                        try:
                            assert(count==1)
                        except:
                            print("Error: Incorrect Shape!")
                            exit(1)
                        self.basisList.append(newBasis)
                        pass
                    elif isEqualVertex(segmentSource.v2, segmentDest.v2):
                        if segmentSource.portIndex==segmentDest.portIndex:
                            newBasis = Basis(segmentSource.v1,\
                                             segmentSource.v2,\
                                             segmentDest.v1,\
                                             segmentSource.portIndex)
                            count+=1
                            pass
                        else:
                            newBasis = Basis(segmentSource.v1,\
                                             segmentSource.v2,\
                                             segmentDest.v1,\
                                             0)
                            count+=1
                            pass
                        segmentSource.contiguous.append(segmentDest.index)
                        segmentDest.contiguous.append(segmentSource.index)
                        try:
                            assert(count==1)
                        except:
                            print("Error: Incorrect Shape!")
                            exit(1)
                        self.basisList.append(newBasis)
                        pass
                    pass
                continue
            continue
        self.export()
        pass
    def export(self):
        try:
            file = open("Mesh/basis.dat", "w")
            for basis in self.basisList:
                file.write("{:21.14E} {:21.14E} {:21.14E} ".format(basis.rm.x, basis.rm.y, basis.rm.z))
                file.write("{:21.14E} {:21.14E} {:21.14E} ".format(basis.rn.x, basis.rn.y, basis.rn.z))
                file.write("{:21.14E} {:21.14E} {:21.14E} ".format(basis.rp.x, basis.rp.y, basis.rp.z))
                file.write("{:3d}\n".format(basis.portIndex))
                continue
            file.close()
            file = open("Mesh/basisInfo.dat", "w")
            file.write(f"{len(self.basisList)}")
            file.close()
        except:
            print("Error: Can't Open Basis File!")
        pass