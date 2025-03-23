import numpy as np
import solver

def FrcTrqArray(bodies,intbodi,intbodj,intbodsi,intbodsj,Fi,Ti):
    #intbodX --> body index for body i and j; 
    #intbodsX --> index for the point at which force is applied
    #Fx --> Force vector on body i and j
    #Tx --> Torque on body i and j
    F =  np.zeros((solver.tnb*3,1))
    if (intbodi == 0):
        Thj =  -bodies[intbodj].xy_global_joints[:,intbodsj][0]*Fi[0,1]+bodies[intbodj].xy_global_joints[:,intbodsj][1]*Fi[0,0]
        F[intbodj*3-1]=-Ti+Thj;F[intbodj*3-2]=-Fi[0,1];F[intbodj*3-3]=-Fi[0,0]    
    elif(intbodj == 0):
        Thi =  bodies[intbodi].xy_global_joints[:,intbodsi][0]*Fi[0,1]-bodies[intbodi].xy_global_joints[:,intbodsi][1]*Fi[0,0]
        F[intbodi*3-1]=Ti+Thi;F[intbodi*3-2]=Fi[0,1];F[intbodi*3-3]=Fi[0,0]
    else: 
        Thi =  bodies[intbodi].xy_global_joints[:,intbodsi][0]*Fi[0,1]-bodies[intbodi].xy_global_joints[:,intbodsi][1]*Fi[0,0]
        Thj =  -bodies[intbodj].xy_global_joints[:,intbodsj][0]*Fi[0,1]+bodies[intbodj].xy_global_joints[:,intbodsj][1]*Fi[0,0]
        F[intbodi*3-1]=Ti+Thi;F[intbodi*3-2]=Fi[0,1];F[intbodi*3-3]=Fi[0,0]
        F[intbodj*3-1]=-Ti+Thj;F[intbodj*3-2]=-Fi[0,1];F[intbodj*3-3]=-Fi[0,0]
    return F

def knmtcs(bodies,intbodi,intbodj,intbodsi,intbodsj,x):
    #distance and velocities between 2 bodies
    #___________________________________________________________________________________________________
    pdist = np.transpose(bodies[intbodi].xy_global_center)+ np.transpose(bodies[intbodi].xy_global_joints[:,intbodsi])-np.transpose(bodies[intbodj].xy_global_center)-np.transpose(bodies[intbodj].xy_global_joints[:,intbodsj])
    if (intbodi == 0):
        xj = x[0,((solver.tnb*3)+(intbodj-1)*3):((solver.tnb*3)+3+(intbodj-1)*3)]
        pvel = -[[xj[0,0],xj[0,1]]]-xj[0,2]*bodies[intbodj].xy_global_jointrot[:,intbodsj]
    elif(intbodj == 0):
        xi = x[0,((solver.tnb*3)+(intbodi-1)*3):((solver.tnb*3)+3+(intbodi-1)*3)]
        pvel = [[xi[0,0],xi[0,1]]]+xi[0,2]*bodies[intbodi].xy_global_jointrot[:,intbodsi]
    else: 
        xi = x[0,((solver.tnb*3)+(intbodi-1)*3):((solver.tnb*3)+3+(intbodi-1)*3)]	
        xj = x[0,((solver.tnb*3)+(intbodj-1)*3):((solver.tnb*3)+3+(intbodj-1)*3)]
        pvel = [[xi[0,0],xi[0,1]]]+xi[0,2]*bodies[intbodi].xy_global_jointrot[:,intbodsi]-[[xj[0,0],xj[0,1]]]-xj[0,2]*bodies[intbodj].xy_global_jointrot[:,intbodsj] 
    return pdist,pvel 
