import numpy as np 

class body_coordinates():
    def __init__(self, index = 0, xy_global_center = [[0.0],[0.0]], xy_local_shape =[[0.0],[0.0]], xy_local_joints = [[0.0],[0.0]], xy_local_unit = [[0.0],[0.0]], angle = 0.0, mass = 0.0, inertia = 0.0): 
        "shape points and local joint points are defined with respect to the center of mass" 
	A90 = [[0.0, -1.0],[1.0, 0.0]]
        self.index =  index
        self.xy_global_center = xy_global_center
	self.angle = angle
        self.xy_local_shape = xy_local_shape
	self.xy_local_joints = xy_local_joints
	self.xy_local_unit = xy_local_unit
	self.mass = [mass,mass,inertia]
	self.xy_jointrot = np.matmul(A90,self.xy_local_joints)
	self.xy_local_unitrot = np.matmul(A90, self.xy_local_unit)

    def BC_trans(self, new_angle, new_position): 
	self.angle = new_angle
	self.xy_global_center = new_position
	self.A = np.array([[np.cos(self.angle), -np.sin(self.angle)],
	                   [np.sin(self.angle), np.cos(self.angle)]]) 
        self.xy_global_shape = self.xy_global_center + np.matmul(self.A,self.xy_local_shape)
	self.xy_global_unitrot = np.matmul(self.A, self.xy_local_unitrot) 
	self.xy_global_joints = np.matmul(self.A,self.xy_local_joints) 
        self.xy_global_jointrot = np.matmul(self.A,self.xy_jointrot) 
        #self.xy_global_shape = np.matmul(self.A,self.xy_local_shape)
class joints(): 
    def __init__(self, joint_type, index, body_i, body_j):
	self.joint_type = joint_type
	self.index = index
	self.body_i = body_i #contains [body_index, joint_on_body_number]
	self.body_j = body_j #contains [body_index, joint_on_body_number]

def BC_constraints(body_i, num_jnt_i, body_j, num_jnt_j, joint_type, qi, qj, tot_num_bod):
    #body_i.BC_trans(qi[2], np.transpose([[qi[0], qi[1]]]))
    #body_j.BC_trans(qj[2], np.transpose([[qj[0], qj[1]]]))

    pt_i = body_i.xy_global_jointrot[:,num_jnt_i] #- np.transpose(body_i.xy_global_center)
    pt_i = np.matrix(pt_i)
    pt_i = np.transpose(pt_i)

    pt_j = body_j.xy_global_jointrot[:,num_jnt_j] #- np.transpose(body_j.xy_global_center)
    pt_j = np.matrix(pt_j)
    pt_j = np.transpose(pt_j)

    I = np.eye(2)

###########################################################################################
				#JOINT TYPE I,J JACOBIANS
###########################################################################################
    if  joint_type == "revolute revolute":
	Jaci = []
        Jacj = []
    elif joint_type == "translational":
        Jaci = []
        Jacj = []
    else:
	"revolute"
	"Constraint Jacobian:"  
	dof_removed = 2       
	Jaci = np.concatenate((-I, -pt_i), axis = 1)
	Jacj = np.concatenate((I, pt_j), axis = 1)
        "Accleration RHS:"
        #print(qi)
        #print(qi[0,5])
        gammai = -qi[0,5]*qi[0,5]*body_i.xy_global_joints[:,num_jnt_i]
	gammaj = -qj[0,5]*qj[0,5]*body_j.xy_global_joints[:,num_jnt_j]
        #gammaij = gammai[:,num_jnt_i] - gammaj[:,num_jnt_j]
	gammaij = gammai - gammaj


###########################################################################################
				#CONCATENATE JACOBIANS#
###########################################################################################
    
    body_diff = body_j.index-body_i.index
    if body_i.index == 0:
        if body_j.index == tot_num_bod:
	    zero1 = np.zeros((dof_removed,(tot_num_bod-1)*3))
	    Jacij = np.concatenate((zero1,Jacj), axis = 1)
        elif body_j.index == 1:
	    zero1 = np.zeros((dof_removed,(tot_num_bod-1)*3))
	    Jacij = np.concatenate((Jacj,zero1), axis = 1) 
        else: 
	    zero1 = np.zeros((dof_removed,(body_j.index-1)*3)) 
	    zero2 = np.zeros((dof_removed,(tot_num_bod-1)*3))
	    Jacij = np.concatenate((zero1,Jacj),axis = 1)
	    Jacij = np.concatenate((Jacij,zero2),axis = 1)
    elif body_j.index == 0:
        if body_i.index == tot_num_bod:
	    zero1 = np.zeros((dof_removed,(tot_num_bod-1)*3))	    
	    Jacij = np.concatenate((zero1,Jaci),axis = 1)
        elif body_i.index == 1:
	    zero1 = np.zeros((dof_removed,(tot_num_bod-1)*3))
	    Jacij = np.concatenate((Jaci,zero1),axis = 1) 
        else: 
	    zero1 = np.zeros((dof_removed,(body_i.index-1)*3)) 
	    zero2 = np.zeros((dof_removed,(tot_num_bod-1)*3))
	    Jacij = np.concatenate((zero1,Jaci))
	    Jacij = np.concatenate((Jacij,zero2)) 
###
    elif (body_diff >= 1):
        zero3 = np.zeros((dof_removed,((tot_num_bod-body_j.index)*3)))
        if body_i.index == 1:
            Jacij = Jaci
        else:
	    zero1 = np.zeros((dof_removed,(body_i.index-1)*3))        
            Jacij = np.concatenate((zero1, Jaci),axis=1)

        if body_diff == 1:  
            zero1 = np.zeros((dof_removed,(body_i.index-1)*3))          
            Jacij = np.concatenate((Jaci,Jacj), axis = 1) 
            Jacij = np.concatenate((zero1,Jacij), axis = 1)  
            Jacij = np.concatenate((Jacij,zero3), axis = 1)
	    #print(Jacij)
        else: 
	    zero2 = np.zeros((dof_removed, body_diff*3)) 
	    Jacij = np.concatenate((Jacij,zero2),axis = 1)
            Jacij = np.concatenate((Jacij,Jacj),axis = 1)

        #if body_j.index != tot_num_bod:
            #zeros3 = np.zeros((dof_removed,(tot_num_bod-body_j.index)*3))
	    #Jacij = np.concatenate((Jacij,zero3), axis = 1)
	    
    elif (body_diff <= -1): 
        if body_j.index == 1:
            Jacij = Jacj
        else:
	    zero1 = np.zeros(dof_removed,(body_j.index-1)*3)        
            Jacij = np.concatenate((zero1, Jacj),axis=1)

        if body_diff == -1:
            Jacij = np.concatenate((Jacij,Jaci),axis = 1)    
        else: 
	    zero2 = np.zeros((dof_removed, body_diff*3)) 
	    Jacij = np.concatenate((Jacij,zero2),axis = 1)
	    Jacij = np.concatenate((Jacij,Jaci),axis = 1)

        if body_j.index != tot_num_bod:
            zeros3 = np.zeros((dof_removed,(tot_num_bod-body_j.index)*3))
	    Jacij = np.concatentate((Jacij,zero3), axis = 1)     

    else:
        Jacij = "Oops, something is wrong!"        
   
    return Jacij, gammaij



     

	           
