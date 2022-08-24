import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
#from IPython.display import IFramepip

def E0(x):
    return np.zeros([1,3])

def B0(x):
    return np.zeros([1,3])

def E1(x):
    return np.array([1.,1.,1.])

def B1(x):
    return np.array([0.,2.,10.])

def Bz(x):
    return np.array([0.,0.,100.])
def push(part, dt):
    part[0] += part[1]*dt
    return part

def Boris(part, E, B, dt,q=1.,m=1.):
    #calculate tt, bisecting angle
    tt= q*B*dt/(2*m)
    #kinetic energy should not change here
    KE = np.linalg.norm(part[1])**2
    #print('initial KE=',KE)
    #calculate s from t
    s= 2*tt / (1+np.linalg.norm(tt)**2)
    #calculate various velocities
    v_minus = part[1] + q/(2*m)*E*dt
    v_prime = v_minus + np.cross(v_minus, tt)
    v_plus = v_minus + np.cross(v_prime, s)
    #v n+1/2
    new_v = v_plus + q/(2*m)*E*dt
    new_part = np.zeros([2,3])
    new_part[0]=part[0]
    new_part[1]=new_v
    KE = np.linalg.norm(new_part[1])**2
    #print('second KE =', KE)
    return new_part
    #return

def simulation(part, Efunc=E0, Bfunc=B0, dt=.01, tot_time=10, q=1., m=1., show_steps=False):
    #t= np.arange(0,tot_time,dt)
    t=[]
    posarr = np.array([])
    velarr = np.array([])
    #push velocity back half step
    E = Efunc(part[0])
    B = Bfunc(part[0])
    #We would like to ensure that u, the magnetic moment, is constant
    #u := W_perp/B = 1/2 * (v_x**2 + v_y**2)
    u = B
    part = Boris(part, E, B, -.5*dt)
    i=0
    while i < tot_time/dt: 
        E = Efunc(part[0])
        B = Bfunc(part[0])
        part = Boris(part,E,B,dt,q,m)
        posarr = np.append(posarr, part[0])
        velarr = np.append(velarr, part[1])
        push(part, dt)
        #code to print positions every so often for debugging
        if show_steps:
            if i>0:
                if i*5*dt/tot_time%1==0:
                    print('timestep ',i*dt," of ", tot_time)
                    #print(part[0])
                    #print('u = ', 1/2*(part[1,0]**2 + part[1,1]**2)/np.linalg.norm(B))
                    KE = np.linalg.norm(part[1])**2
                    print('KE=',KE)
        i+=1
        t.append(i*dt)
        if sum(part[0]**2)**.5>20:
            print('simulation broke at time ',i*dt)
            print('final position was', part[0], ' which is ',sum(part[0]**2)**.5, 'away from origin')
            break
    posarr = np.reshape(posarr, [len(t),3])
    velarr = np.reshape(velarr, [len(t),3])
    print('initial', posarr[0])
    print('next few', posarr[30])
    #return positions and velocities
    return (posarr, velarr)



#I need to fix this with the updated s
#fixed!!!
def func(num, dataSet, line):
    #x and y axis, i.e, the first and second arrays in the first dim, until num
    #old version
    #line.set_data(dataSet[0:2, :num])
    #new
    space = 10
    if num>space:
 #   if False:
        
        line.set_data((dataSet[num-space:num, 0],dataSet[num-space:num, 1]))
        line.set_3d_properties(dataSet[num-space:num,2])
    else:
        line.set_data((dataSet[:num, 0],dataSet[:num, 1]))
        line.set_3d_properties(dataSet[:num,2])
    # line.set_data((dataSet[:num, 0],dataSet[:num, 1]))
    # line.set_3d_properties(dataSet[:num,2])
    return (line)
def B_bottle(pos, mu_0=100, dipole_distance=10, mu = np.array([0,0,5])):
    x,y,z = pos;
    #set the position of the top dipole r prime
    rp_top = (0.0, 0.0, dipole_distance)
    #calculate r vector
    r_top = np.array([x,y,z]) - rp_top
    #print(r_top)
    #find magnitude of r
    rmag = sum(r_top**2)**.5
    #print(rmag)
    #B(r)
    Btop = mu_0/(4*np.pi)*(3*r_top * np.dot(mu, r_top)/np.power(rmag,5) - mu/np.power(rmag,3))
    #repeat for bottom dipole
    rp_bot = (0.0, 0.0, -dipole_distance)
    #calculate r vector
    r_bot = np.array([x,y,z]) - rp_bot
    #find magnitude of r
    rmag = sum(r_bot**2)**.5
    #B(r)
    Bbot = mu_0/(4*np.pi)*(3*r_bot * np.dot(mu, r_bot)/rmag**5 - mu/rmag**3)
    return (Btop + Bbot)

def dipole(pos, dpos=np.array([0,0,10]), m=np.array([0,0,1])):
    r = pos-dpos
    return (np.dot(m,r)*r-m)/np.linalg.norm(r)**3
def B_d(pos, dd=10, B=60):
    dpos = np.array([0,0,dd])
    Btop = dipole(pos, dpos)
    Bbot = dipole(pos, -dpos)
    return B*(Btop + Bbot)

def coil(pos,x_p,R):
    x=pos-x_p
    return .5* R**2 * np.power((R**2+x**2),-3/2)

def B_h(pos, dd=10., R=8):
    Btop = coil(pos, np.array([0,0,dd]),R)
    Bbot = coil(pos, np.array([0,0,-dd]),R)
    return Btop + Bbot

def create_plot(m=1.,q=1.,dt=.001,tot_time=10., part=np.zeros([2,3]),E=E1, B=B1 ,s=10, filename='anim', l=20,steps=False):
    print('start is', part[0])
    (pos, vel) = simulation(part, E, B, dt, tot_time,q=q,m=m, show_steps=steps)
    print('big size is ',len(pos))
    # s = int(len(pos)/N)
    # print(s)
    # big_dataSet = np.reshape(pos, [3,len(pos)])
    # if s>2:
    #     dataSet = big_dataSet[:, :-1:s]
    # else:
    #     dataSet = big_dataSet[:, :-1:2]
    #numDataPoints = int(len(dataSet[0]))
    
    s_pos = pos[::s]
    print('ndp=',len(s_pos))
    numDataPoints = int(len(s_pos))
    
    fig = plt.figure()
    ax = Axes3D(fig)
    # NOTE: Can't pass empty arrays into 3d version of plot()
    line = plt.plot(s_pos[:,0], s_pos[:,1], s_pos[:,2], lw=2, c='g')[0] # For line plot
    ax.scatter(0,0,10)
    ax.scatter(0,0,-10)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Mytitle')
    
    ax.set_xlim([-l,l])
    ax.set_ylim([-l,l])
    ax.set_zlim([-l,l])

    # Creating the Animation object
   # line_ani = animation.FuncAnimation(fig, func, frames=numDataPoints, fargs=(s_pos,line), interval=50, blit=False)
    line_ani = animation.FuncAnimation(fig, func, frames=numDataPoints, fargs=(s_pos,line), interval=50, blit=False)

    writergif = animation.PillowWriter(fps=30)
    s_filename = filename + '.gif'
    line_ani.save(s_filename,writer=writergif)

    # resized output IFrame
    #IFrame(src=s_filename, width=600, height=450)
    print('done')
    return s_pos

test =np.zeros([2,3])
test[0,0]=5.
test[0,2]=-3.
test[1,1]=-10.1
test[1,2]=.1
test[1,0]=.1
d = create_plot(part = test,E=E0, B=B_d, s=10, tot_time=20.,l=10, dt=.004,steps=True, filename='new_animation')
print('gif saved!')