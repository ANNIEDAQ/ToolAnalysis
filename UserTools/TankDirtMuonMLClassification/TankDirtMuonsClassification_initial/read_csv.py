import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as scp
import numpy as np

df = pd.read_csv ('tank_muons_data_1000evts.csv')
print (df.head())
#print("All columns are: ", df.columns.values.tolist())
print("")
#print(df[:10])
#df.to_csv("test.csv",index=False) #write dataframe to csv file

 #---- My Plots:
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
fig = plt.figure()
ax3 = fig.add_subplot(111, projection='3d')
#ax.plot_surface(x_coord, y_coord, z_coord, cmap=cm.coolwarm, linewidth=0, antialiased=False);
xval = df.loc[0, "pmtx_0":"pmtx_5999"]
yval = df.loc[0, "pmty_0":"pmty_5999"]
zval = df.loc[0, "pmtz_0":"pmtz_5999"]
#ax3.scatter(xs=xval, ys=yval,zs=zval, zdir='z', s=qval*2, c='m');
ax3.scatter(xs=xval, ys=yval,zs=zval,c='m');
ax3.set_xlabel("x [mm]")
ax3.set_ylabel("y [mm]")
ax3.set_zlabel("z [mm]")
plt.show()

fig, ax0 = plt.subplots()
ax0 = fig.add_subplot(111)
xval = df.loc[0, "pmtx_0":"pmtx_5999"]
yval = df.loc[0, "pmty_0":"pmty_5999"]
ax0.scatter(xval, yval);
#df.plot(x=df.loc[0, "pmtx_0":"pmtx_5999"] ,y=df.loc[0, "pmty_0":"pmty_5999"], ax=ax0, kind='scatter', c='m',grid=True)
##ax0.set_xlim(-3700,3700)
##ax0.set_ylim(-3700,3700)
plt.show()
