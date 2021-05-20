import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
filelist=[]

#liimited_filelist = []
#iC_square=np.loadtxt('IC.txt')

#d1 = np.loadtxt('sol_10.txt')
x = np.linspace(0, 1, 100)
j = 1
for i in range(0,199):
    #if i % 10 == 0:
    filelist.append("sol_{}.txt".format(i))
    #liimited_filelist.append("FROMM_UNLIMITED_{}.txt".format(i))

#colour=iter(cm.rainbow(np.linspace(0,20,860)))
#plt.plot(x, d1, marker = 'o', markersize = 1, linewidth = 2, color = 'blue', label = 't = 0.0 s')
#plt.plot(x, iC_square, marker = 'o', markersize = 1.6, linewidth = 2, color = 'darkblue', label = 'CFL = 0.5')
colour=iter(cm.rainbow(np.linspace(0, 20,100)))

#IC = np.loadtxt("IC.txt")
#plt.plot(x, IC, linewidth = 1, color = "darkblue", label = 'CFL = 0.5')
for fname  in filelist:
    
    #c = next(colour)
    #c=next(colour)
    plt.figure(figsize=(8,8))
    data=np.loadtxt(fname)
   
    plt.plot(x, data, linewidth = 1, color = 'blue', label = 'CFL = 0.9')

    plt.ylabel('u (m/s)', fontsize=15)
    plt.xlabel('h (m)', fontsize=15)
    plt.ylim(0, 1.2)
    plt.xlim(0, 1)
    plt.title("Gaussain IC 1D Lax Wendroff No Limited")
    plt.legend(loc = 'upper right')
    plt.savefig('ShallowWater1{}.png'.format(j))
    plt.clf()
    plt.close()

    j += 1

#plt.legend()
plt.show()

'''
import numpy as np
import matplotlib.pyplot as plt

filelist=[]

liimited_filelist = []


x = np.linspace(-1, 1, 200)
j = 1
for i in range(1,100):
    filelist.append("FVM_1S_B&amp;W_{}.txt".format(i))
    liimited_filelist.append("BeamW_Limited_{}.txt".format(i))

for fname, fname1 in zip(filelist, liimited_filelist):
    plt.figure(figsize=(12,6))
    data=np.loadtxt(fname)
    data1=np.loadtxt(fname1)
    #initial_data = np.loadtxt(initial_file[0])
    #initial_data = np.loadtxt(filelist[0])
    plt.plot(x, data, color = 'orange', linewidth = 0.7, label = 'B&W unlimited')
    plt.plot(x, data1, marker = 'o', markersize = 4, linewidth = 1.3, color = 'blue', label = 'B&W Method 2nd Order')
    #plt.plot(x, data, marker = 'o', markersize = 4, linewidth = 2, color = 'darkorange', label = 'B&W method Van\nAlbada Slope Limited')
    plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
    plt.title('1D Linear Convection 2nd Order Beam & Warming FVM Analysis ', fontsize=15)
    plt.xlabel('u (m/s)', fontsize=15)
    plt.ylabel('time t(s)', fontsize=15)
    plt.ylim(-0.4, 1.4)
    plt.xlim(-1.1, 1.1)
    #plt.title("1D Linear Convection FVM")
    plt.legend(loc = 'upper right')
    plt.savefig('B&W{}.png'.format(j))
    plt.clf()
    plt.close()
    j += 1

#plt.legend()
plt.show()




'''




# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# filelist=[]

# #liimited_filelist = []
# #iC_square=np.loadtxt('IC.txt')

# #d1 = np.loadtxt('sol_10.txt')
# x = np.linspace(0, 1, 100)
# j = 1
# for i in range(0, 39):
#     if i %  6 == 0:
#         filelist.append("sol_{}.txt".format(i))
#     #liimited_filelist.append("FROMM_UNLIMITED_{}.txt".format(i))

# #colour=iter(cm.rainbow(np.linspace(0,20,860)))
# #plt.plot(x, d1, marker = 'o', markersize = 1, linewidth = 2, color = 'blue', label = 't = 0.0 s')
# #plt.plot(x, Solution_0, marker = 'o', markersize = 1.6, linewidth = 2, color = 'darkblue', label = 'CFL = 0.5')
# colour=iter(cm.rainbow(np.linspace(0, 9,58)))
# plt.figure(figsize=(8,8))
# #IC = np.loadtxt("IC.txt")
# #plt.plot(x, IC, linewidth = 1, color = "darkblue", label = 'CFL = 0.5')
# #plt.plot(-2, 5, marker = 'o', markersize = 2, linewidth = 1, label = "CFL = 0.9")
# for fname  in filelist:
    
#     #plt.figure(figsize=(8,8))
#     c = next(colour)
#     #c=next(colour)
#     data=np.loadtxt(fname)
#     #data1=np.loadtxt(fname1)
#     #initial_data = np.loadtxt(initial_file[0])
#     #initial_data = np.loadtxt(filelist[0])
#     plt.plot(x, data, linewidth = 1.5, color = c)

#     #plt.plot(x1, initial_data, color = 'orange', linewidth = 0.7, label = 'Initial Condition')
#     #plt.plot(x, data, marker = 'o', markersize = 2.1, linewidth = 1.3, color = 'blue', label = 'FROMM Method 2nd Order')
#     #plt.plot(x, data, marker = 'o', markersize = 4, linewidth = 2, color = 'darkorange', label = 'Fromm method Van\nAlbada Slope Limited')
#     #plt.plot(-1,1, color = 'black', label = 'timestep {}'.format(j))
#     #plt.title('1D Inviscid Burgers Sine wave ICs ', fontsize=15)
#     plt.xlabel('x (m)', fontsize=15)
#     plt.ylabel('h (m)', fontsize=15)
#     plt.ylim(0, 1.2)
#     plt.xlim(0, 1)
#     plt.title("Lax Wendroff Van Albada, CFL = 0.9, Time = 0.2 s".format(j))
#     plt.legend(loc = 'upper right')
#     plt.savefig('LWVanAlbada_LimiterH{}.png'.format(j))
#     #plt.clf()
#     #plt.close()


#     j += 1

# #plt.legend()
# plt.show()


