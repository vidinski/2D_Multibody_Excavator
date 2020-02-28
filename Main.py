import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import scipy.integrate as integrate
import solver
from body_coordinates import body_coordinates, BC_constraints, joints

PI = np.pi; cos45 = np.cos(PI/4.)
solver.bodies = []
go = []
solver.tnb = 3
dt = 0.1

tspan = [0., 10.0]
tmr = np.arange(0.0, tspan[1], dt)

###########################################################################################
					#SETUP BODIES
###########################################################################################

####GROUND#####
body0 = body_coordinates(0, #index
			 np.transpose([[0.0, 0.0]]), #center
			 np.transpose([[0.0, 0.0]]), #shape
			 np.transpose([[1.2, 0.0],[0.0,0.6]]), #joints
			 np.transpose([[0.0, 1.0], [1.0,0.0]]), #unit vectors
			 0.0, #angle
			 100000.0, #mass
			 100000.0) #inertia 
solver.bodies.append(body0)
solver.ground0true = 1		
##_______________________________________________________________________________________##
th1 = 0.333*PI
body1 = body_coordinates(1, #index
			 np.transpose([[1.2*np.cos(th1),1.2*np.sin(th1)+1.2]]), #center
			 np.transpose([[1.2,0.0],[1.2+0.3*np.cos(PI/4.), 0.3-0.3*np.sin(PI/4.)],
			              	[1.5, 0.3],[1.2+0.3*np.cos(PI/4.),0.3+0.3*np.sin(PI/4.)],
					[1.2,0.6],[-1.2,0.6],
					[-1.2-0.3*np.cos(PI/4.),0.3+0.3*np.sin(PI/4.)],[-1.5,0.3],
					[-1.2-0.3*np.cos(PI/4.),0.3-0.3*np.sin(PI/4.)],[-1.2,0.0]]), #shape
			 np.transpose([[-1.2,0.0],[1.2,0.0],[0.8,-0.3],[-0.8,0.3]]), #joints
			 np.transpose([[0.0,1.0],[1.0, 0.0]]), #unit vectors
			 th1, #angle
			 300.0, #mass
			 300.0/12.0*(0.6**2+3.6**2)) #inertia
solver.bodies.append(body1)
go.append(patches.Polygon(np.array(np.transpose(solver.bodies[1].xy_local_shape)),edgecolor='black', facecolor='orange'))

##_______________________________________________________________________________________##
th2 = 0.5*PI; th12 = th1-th2
body2 = body_coordinates(2, #index
			 np.transpose([[0.6*np.cos(th12),0.6*np.sin(th12)]])+np.transpose([[2.4*np.cos(th1),2.4*np.sin(th1)+1.2]]), #center
			 np.transpose([[1.5,0.0],[1.5+0.3*np.cos(PI/4.), 0.3-0.3*np.sin(PI/4.)],
			              	[1.8, 0.3],[1.5+0.3*np.cos(PI/4.),0.3+0.3*np.sin(PI/4.)],
					[1.5,0.6],[-1.5,0.6],
					[-1.5-0.3*np.cos(PI/4.),0.3+0.3*np.sin(PI/4.)],[-1.8,0.3],
					[-1.5-0.3*np.cos(PI/4.),0.3-0.3*np.sin(PI/4.)],[-1.5,0.0]]), #shape
			 np.transpose([[-0.6,0.0],[1.5,0.0],[-1.0,0.0],[-0.8,0.3]]), #joints
			 np.transpose([[0.0,1.0],[1.0, 0.0]]), #unit vectors
			 th12,   #angle
			 400.,      #mass
			 400.0/12.0*(0.6**2+3.6**2)) #inertia
solver.bodies.append(body2)
go.append(patches.Polygon(np.array(np.transpose(solver.bodies[2].xy_local_shape)),edgecolor='black', facecolor='orange'))


##_______________________________________________________________________________________##
body3 = body_coordinates(3, #index
			 body2.xy_global_center+np.transpose([[1.5*np.cos(th12),1.5*np.sin(th12)]])+np.transpose([[0.3,-0.6]]), #center
			 np.transpose([[0.6,-0.6],[0.6,0.6],[-0.3,0.6],[-0.6,0.0],[-0.3,-0.6],[-0.3,-1.2]]), #shape
			 np.transpose([[-0.3,0.6],[0.6,0.6]]), #joints
			 np.transpose([[0.0,1.2],[1.0, 0.0]]), #unit vectors
			 0.0,   #angle
			 100.,      #mass
			 100./12.0) #inertia
solver.bodies.append(body3)
go.append(patches.Polygon(np.array(np.transpose(solver.bodies[3].xy_local_shape)),edgecolor='orange', facecolor='black'))




###########################################################################################
					#SETUP ANIMATION WINDOW
###########################################################################################
plt.figure(2)
fig = plt.figure(1)
ax1 = fig.add_subplot(111)

time_template = 'time = %.1fs'
time_text = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)

ax1.set_ylim(-2,8)
ax1.set_xlim(-4,6)
background1 = patches.Rectangle([-50., 0.], 100.0*2., 50.0*2., 0.0, color = 'blue', alpha = 0.35)
background2 = patches.Rectangle([-50., 0.], 100.0*2., -50.0*2., 0.0, color = 'green', alpha = 0.5)
background3 = patches.Polygon(np.array([[0.,0.9],[0.0,2.4],[-0.6, 3.0],[-1.2,3.0],[-1.8,2.1],[-2.4,1.5],[-2.4,0.9]]), edgecolor='black', facecolor='orange')
background4 = patches.Polygon(np.array([[0.9,0.0],[0.9+0.3*np.cos(PI/4.),0.3-0.3*np.sin(PI/4.)],
			                [1.2, 0.3],[0.9+0.3*np.cos(PI/4.),0.3+0.3*np.sin(PI/4.)],
					[0.9,0.6],[-3.0,0.6],
					[-3.0-0.3*np.cos(PI/4.),0.3+0.3*np.sin(PI/4.)],[-3.3,0.3],
					[-3.0-0.3*np.cos(PI/4.),0.3-0.3*np.sin(PI/4.)],[-3.,0.0]]), color='black')
background5 = patches.Rectangle([-1.5, 0.6], 0.9, 0.3, 0.0, edgecolor = 'orange', facecolor = 'black')

background1 = ax1.add_patch(background1)
background2 = ax1.add_patch(background2)
background3 = ax1.add_patch(background3)
background4 = ax1.add_patch(background4)
background5 = ax1.add_patch(background5)

for j in range (1,len(solver.bodies)):
    jj = j-1
    solver.bodies[j].BC_trans(solver.bodies[j].angle,solver.bodies[j].xy_global_center)
    go[jj].set_xy(np.transpose(solver.bodies[j].xy_global_shape))
    patch = ax1.add_patch(go[jj])

###########################################################################################
					#SETUP JOINTS
###########################################################################################
solver.joint_list = []

##_______________________________________________________________________________________##
joint1 = joints("revolute", 1, [body0.index, 1], [body1.index,0])

solver.joint_list.append(joint1)

##_______________________________________________________________________________________##
joint2 = joints("revolute", 2, [body1.index, 1], [body2.index,0])

solver.joint_list.append(joint2)

##_______________________________________________________________________________________##
joint3 = joints("revolute", 3, [body2.index, 1], [body3.index,0])

solver.joint_list.append(joint3)


###########################################################################################
			#SETUP EQUATIONS OF MOTION AND SOLVE 
###########################################################################################
M = solver.bodies[1].mass
x0 = np.concatenate((solver.bodies[1].xy_global_center, np.matrix(solver.bodies[1].angle)), axis = 0)
if len(solver.bodies) > 2:
    for i in range (2,len(solver.bodies)):
        M = np.concatenate((M, solver.bodies[i].mass), axis = 0)
        x0 = np.concatenate((x0,np.concatenate((solver.bodies[i].xy_global_center, np.matrix(solver.bodies[i].angle)), axis = 0) ), axis = 0)


M = np.diag(M)
solver.invM = np.linalg.inv(M)
x0 = np.concatenate((np.transpose(x0),np.zeros([1,3*solver.tnb])), axis = 1)

###ADD STATES###
x0 = np.concatenate((x0,np.zeros([1,12])), axis = 1)
x0 = np.array(x0) 

x0[0,20] = 150e+05; x0[0,21] = 150e+05 

x0[0,24] = 150e+05; x0[0,25] = 150e+05

x0[0,28] = 150e+05; x0[0,29] = 150e+05

solver.leg = np.transpose(solver.bodies[1].xy_global_center) + np.transpose(solver.bodies[1].xy_global_joints[:,0])
solver.anchor = np.matrix([[0.,0.]]) 
solver.flag1 = 1
solver.flag2 = 0
x = integrate.solve_ivp(solver.solveSys, tspan, x0[0], method='RK45',t_eval = tmr)
###########################################################################################
					#ANIMATION
###########################################################################################

def init():
    return [] #patch

def animate(i):
    time_text.set_text(time_template % (i*dt))
    q = np.zeros([solver.tnb,6])   
    x_anim = x.y[:,i] 
    x_anim = np.matrix(x_anim)  
    #print(x_anim)
    for j in range (1,len(solver.bodies)):
        jj = j-1
        q[jj,:] = np.concatenate((x_anim[0,3*jj:3*jj+3],x_anim[0,3*(solver.tnb+jj):3*(solver.tnb+jj)+3]),axis=1)
        solver.bodies[j].BC_trans(q[jj,2], np.transpose([[q[jj,0], q[jj,1]]]))
        go[jj].set_xy(np.transpose(solver.bodies[j].xy_global_shape))
        patch = ax1.add_patch(go[jj]) 
    return [] #patch,


ani = animation.FuncAnimation(fig, animate, np.arange(1, len(np.transpose(x.y))), interval=1, blit=True, init_func=init)

plt.figure(2)
plt.plot(x.t,np.array(x.y[18,:]))
plt.plot(x.t,np.array(x.y[22,:]))
plt.plot(x.t,np.array(x.y[26,:]))
plt.title('Valve Position vs Time')
plt.xlabel('time (s)')
plt.ylabel('Valve Position (m)')
plt.grid(True)
plt.legend(('vcmd1','vcmd2','vcmd3'))
plt.show()












