import numpy as np

# %% Return the Element with Adequate Nodes

def getElement(nbrNode):
    match nbrNode:

        case 2: return Line()
        case 3: return Triangle()
        case 4: return Quadrangle()

# %% Linear Line Finite Element

class Line(object):

    def evaluate(self,parm):
        return np.array([(1-parm[0])/2,(1+parm[0])/2])

    # Position of a parametric point in the reference space
    
    def getPosition(self,node,parm):
        return np.squeeze(self.evaluate(parm).dot(node))
    
    # Projection of a point in the reference space

    def projection(self,node,pos):
    
        A = np.diff(node,axis=0)/2
        B = np.array(pos-np.sum(node,axis=0)/2)
        return np.linalg.lstsq(np.transpose(A),B,-1)[0]

    # Distance between a point and the projection

    def distance(self,parm,node,pos):

        if abs(parm)>1.001: return np.inf
        return np.linalg.norm(self.getPosition(node,parm)-pos)

# %% Linear Triangle Finite Element

class Triangle(Line):

    def evaluate(self,parm):
        return np.array([1-parm[0]-parm[1],parm[0],parm[1]])
    
    # Projection of a point in the reference space

    def projection(self,node,pos):

        B = np.array(pos-node[0])
        A = [node[1]-node[0],node[2]-node[0]]
        return np.linalg.lstsq(np.transpose(A),B,-1)[0]

    # Distance between a point and the projection

    def distance(self,parm,node,pos):

        if any(-0.001>parm) or sum(parm)>1.001: return np.inf
        return np.linalg.norm(self.getPosition(node,parm)-pos)

# %% Linear Quadrangle Finite Element

class Quadrangle(Line):

    def evaluate(self,parm):

        return np.array([
        (1-parm[0])*(1-parm[1])/4,(1+parm[0])*(1-parm[1])/4,
        (1+parm[0])*(1+parm[1])/4,(1-parm[0])*(1+parm[1])/4])
    
    # Gradient of the element shape functions
    
    def grad(self,parm):

        return np.array([
        [(parm[1]-1)/4,(1-parm[1])/4,(1+parm[1])/4,(-parm[1]-1)/4],
        [(parm[0]-1)/4,(-parm[0]-1)/4,(1+parm[0])/4,(1-parm[0])/4]])
    
    # Distance between a point and the projection

    def distance(self,parm,node,pos):

        if any(abs(parm)>1.001): return np.inf
        return np.linalg.norm(self.getPosition(node,parm)-pos)

    # Projection of a point in the reference space

    def projection(self,node,pos):

        residual = np.inf
        parm = np.zeros(2)

        # Newton iterations for parametric coordinates

        while np.any(abs(residual)>1e-12):

            B = self.getPosition(node,parm)-pos
            A = np.atleast_2d(self.grad(parm).dot(node))
            residual = np.linalg.lstsq(np.transpose(A),B,-1)[0]
            parm = parm-residual

        return parm
