import numpy as np
#from scipy.integrate import odeint
from body_coordinates import body_coordinates, BC_constraints, joints
import frconbod

global bodies
global invM
global K
global tnb
global joint_list
global ground0true



######################################################
beta = 1.376*1e+9 	#bulk modulus
rho = 857.0 		#density
Dh = 0.01        	#hydraulic diameter
Ar = 0.1 		#Retract area
Ae = 0.1		#Extend Area
Cd = 0.7		#discharge coefficient
Ps = 300.0*1e+05	#Source Pressure
zetas = 0.707       	#damping spool N/m/s**2
wns = 30.0*2.0*3.14	#Natural frequency valve
ka = 1.77765e+03	#current gain to spool
rl = 1.75		#rod length
Voe = 5e-05
Vor = 1.75*Ae
bl = 1.75
cmd1off = 0.0
g = 9.81
PI = np.pi

def solveSys(t,x):
    #_______________________________________________________________________
    x = np.matrix(x)
    q = np.zeros([tnb,6]); 
    bodies[0].BC_trans(0.0, np.transpose([0.0, 0.0])) 
    D = []
    for i in range (1,len(bodies)):
        j = i-1
        q[j,:] = np.concatenate((x[0,3*j:3*j+3],x[0,3*(tnb+j):3*(tnb+j)+3]), axis = 1)
        bodies[i].BC_trans(q[j,2], np.transpose([[q[j,0], q[j,1]]]))
    #Constraints Jacobian: 
    #_______________________________________________________________________
    if (ground0true == 1):
        [Jacij, gammaij] = BC_constraints(bodies[joint_list[0].body_i[0]],
		          	      joint_list[0].body_i[1],
		   		      bodies[joint_list[0].body_j[0]], 
		   		      joint_list[0].body_j[1],
		   		      joint_list[0].joint_type, 
				      np.matrix([0,0,0,0,0,0]),
                                      np.matrix(q[0,:]), 
    				      tnb)
        D = Jacij
        gamma = np.transpose(np.matrix(gammaij))

        if len(joint_list) > 1: 
           for i in range (1,len(joint_list)):  
	        #print(joint_list[i].body_i[1])
                 [Jacij, gammaij] = BC_constraints(bodies[joint_list[i].body_i[0]],
		          	              joint_list[i].body_i[1],
		   		              bodies[joint_list[i].body_j[0]], 
		   		              joint_list[i].body_j[1],
		   		              joint_list[i].joint_type, 
		   		              np.matrix(q[joint_list[i].body_i[0]-1,:]), 
                                              np.matrix(q[joint_list[i].body_j[0]-1,:]),
                                              tnb)
	         D = np.concatenate((D,Jacij),axis = 0)
	         gammaij = np.transpose(np.matrix(gammaij))
	         gamma = np.concatenate((gamma,gammaij),axis = 0)
    else:
        [Jacij, gammaij] = BC_constraints(bodies[joint_list[0].body_i[0]],
		          	      joint_list[0].body_i[1],
		   		      bodies[joint_list[0].body_j[0]], 
		   		      joint_list[0].body_j[1],
		   		      joint_list[0].joint_type, 
		   		      np.matrix(q[joint_list[0].body_i[0]-1,:]), 
                                      np.matrix(q[joint_list[0].body_j[0]-1,:]),
    				      tnb)
        

        D = Jacij
	gamma = np.transpose(np.matrix(gammaij))

        if len(joint_list) > 1: 
           for i in range (1,len(joint_list)):  
	        #print(joint_list[i].body_i[1])
                 [Jacij, gammaij] = BC_constraints(bodies[joint_list[i].body_i[0]],
		          	              joint_list[i].body_i[1],
		   		              bodies[joint_list[i].body_j[0]], 
		   		              joint_list[i].body_j[1],
		   		              joint_list[i].joint_type, 
		   		              np.matrix(q[joint_list[i].body_i[0]-1,:]), 
                                              np.matrix(q[joint_list[i].body_j[0]-1,:]),
                                              tnb)
	         D = np.concatenate((D,Jacij),axis = 0)
	         gammaij = np.transpose(np.matrix(gammaij))
	         gamma = np.concatenate((gamma,gammaij),axis = 0)
    
    #################   COMPUTE FORCES   ########################################################   

    [F,xd_hyd]= forces(x,t,bodies)

    DT = np.transpose(D)
    DMinv = np.matmul(D,invM)
    DMinvDT = np.matmul(DMinv,DT)
    DMinvF = np.matmul(DMinv, F)
    gamma_F = gamma-DMinvF
    lagr_lambda = np.linalg.solve(DMinvDT,gamma_F)
    qdoubledot = np.matmul(DT,lagr_lambda)
    qdoubledot = qdoubledot+F
    qdoubledot = np.matmul(invM, qdoubledot)   
    
    #extract velocities: 
    xyth_vel = x[0,tnb*3:tnb*6]
    qdoubledot = np.transpose(qdoubledot)
    xdot = np.concatenate((xyth_vel,qdoubledot),axis = 1 )
    
    ############### ADD HYDRAULIC STATES ################## 
    #######################################################
    xdot = np.concatenate((xdot,xd_hyd),axis = 1)
    #######################################################    

    xdot = np.array(xdot)
    return xdot[0]

def forces(x,t,bodies):
    kground = 1e+7
    cground = 1e+5
    g = -9.81*np.matrix([0.0,1.0,0.0])
    anchor1 = [[0.0, 0.0]]
    Fg = g
    for i in range (2,len(bodies)):
        Fg = np.concatenate((Fg,g), axis = 1) 
    Fh = np.zeros((tnb*3,1))
    ### system states ###
    th = np.zeros(((len(bodies)-1)*2,1))
    thd = np.zeros(((len(bodies)-1)*2,1))
    #Step input
    if t <= 1.0 or t>6.0: 
       cmd1 = 0.0    
    else: 
       cmd1 = -0.03 
    
    #################   CYLINDER1   #############################################################
    [distLift,velLift] = frconbod.knmtcs(bodies,1,0,2,1,x)
    pnum = 1 #piston number
    [xd_hyd1,F_hyd1,xp1,xdp1]= forces_hydraulic(x,t,bodies, distLift, velLift, cmd1,pnum)
    Fh1 = (F_hyd1/xp1)*distLift;  
    Fh = frconbod.FrcTrqArray(bodies,1,0,2,1,Fh1,0.0)
    
    #################   CYLINDER2   #############################################################
    [distBoom,velBoom] = frconbod.knmtcs(bodies,2,1,2,3,x)
    cmd2 = 0.75*cmd1
    pnum = 2 #piston number
    [xd_hyd2,F_hyd2,xp2,xdp2]= forces_hydraulic(x,t,bodies, distBoom, velBoom, cmd2,pnum)
    Fh2 = (F_hyd2/xp2)*distBoom 
    FBoom = frconbod.FrcTrqArray(bodies,2,1,2,3,Fh2,0.0)
    Fh = FBoom + Fh
    #################   CYLINDER3   #############################################################
    [distYoke,velYoke] = frconbod.knmtcs(bodies,3,2,1,3,x)
    cmd3 = 1.5*cmd1
    pnum = 3 #piston number
    [xd_hyd3,F_hyd3,xp3,xdp3]= forces_hydraulic(x,t,bodies, distYoke, velYoke, cmd3,pnum)
    Fh3 = (F_hyd3/xp3)*distYoke;  
    #FYokeT = np.matmul(np.transpose(bodies[3].xy_global_jointrot[:,1]),np.transpose(Fh3))
    #Fh[6,0]=Fh3[0,0]; Fh[7,0]=Fh3[0,1]; Fh[8,0]=FYokeT
    #Fh[3,0]=Fh2[0,0]-Fh3[0,0]; Fh[4,0]=Fh2[0,1]-Fh3[0,1];
    Fyoke = frconbod.FrcTrqArray(bodies,3,2,1,3,Fh3,0.0)
    Fh = Fh + Fyoke
    ################# CONCATENATE HYDRAULIC STATES ###############################################
    xd_hydt = np.concatenate((xd_hyd1,xd_hyd2),axis=0)
    xd_hydt = np.concatenate((xd_hydt,xd_hyd3),axis=0)
    xd_hyd = np.transpose(xd_hydt)
    F = np.transpose(Fg)+Fh
    return F,xd_hyd
  
def forces_hydraulic(x,t,bodies,dist,vel,vcmd,pnum):
    #compute x dot for hydraulics 
    xp = np.linalg.norm(dist); 
    xp = np.matmul(dist,1.0/xp*np.transpose(dist)); 
    xdp = np.matmul(vel,1.0/xp*np.transpose(dist)); 
    x1 = x[0,(tnb*6)+4*(pnum-1)]
    x2 = x[0,(tnb*6+4*(pnum-1)+1)]
    x3 = x[0,(tnb*6+4*(pnum-1)+2)] 
    x4 = x[0,(tnb*6+4*(pnum-1)+3)]
    xd_hyd = np.zeros((4,1))
    xd_hyd[0,0] = x2
    xd_hyd[1,0] = ka*vcmd-2.0*zetas*wns*x2-wns**2*x1
    if vcmd >= 0.0:
       xd_hyd[2,0] = beta/(Ae*(xp-rl)+Voe)*np.sign(Ps-x3)*(Cd*Dh*x1*np.sqrt(2.0/rho*np.abs((Ps-x3)))-Ae*xdp)
       xd_hyd[3,0] = beta/(Vor-Ar*(xp-rl))*(-Cd*Dh*x1*np.sqrt(2.0/rho*x4)+Ae*xdp)
    else: 
       xd_hyd[2,0] = beta/(-Ae*(bl-xp+rl)+Voe)*(-Cd*Dh*x1*np.sqrt(2.0/rho*x3)+Ae*xdp)
       xd_hyd[3,0] = beta/(Vor+Ar*(bl-xp+rl))*(Cd*Dh*x1*np.sign(Ps-x4)*np.sqrt(2.0/rho*np.abs(Ps-x4))-Ae*xdp)    

    F_hyd = x3*Ae-Ar*x4; 
    return xd_hyd,F_hyd,xp,xdp


