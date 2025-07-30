#!/usr/bin/python
# Active Selection of Molecules (Short task 01)
# Marcos Quiles / Marinalva Soares

# INPUT1 - XYZ files (folder)
# INPUT2 - Text file - Additional information (Energy_T, Delta_E, etc.)
# INPUT3 - #number of outputs
import os
import time
start=time.time()
import sys
# print(len(sys.argv))
import matplotlib.cm as cm
import matplotlib.pyplot as plt
if len(sys.argv) != 4 and len(sys.argv) != 5:
    print("\nUsage: \n")
    print("\tOption 1: Only structural information (xyz) \n")
    print("\t\t$ python script1.py 1 folder_xyz_files #samples\n\n");
    print("\tOption 2: Structural information (xyz) + Energy_T\n")
    print("\t\t$ python script1.py 2 folder_xyz_files #samples extra_data.txt \n\n");
    print("\tOption 3: Structural information (xyz) + {Energy_T, Delta_E, m_T, ECN, d_av}\n")
    print("\t\t$ python script1.py 3 folder_xyz_files #samples extra_data.txt \n\n");
    exit();


opt = int(sys.argv[1])
if opt not in range(1,4):
    print("Wrong option\n\n");
    exit() 

params = {'eigen': 9999999,
            'n_clusters': int(sys.argv[3])}

import pandas as pd
from glob import glob
from os import makedirs, getcwd, path
from tools import *

baseFolder = glob(str(sys.argv[2])+'/*.xyz')
# print(type(baseFolder))
baseFolder.sort()
dataxyz = []
ids = []
num_files = 1

for fin in baseFolder:
    # print(fin)
    natoms, atomtypes, coords = xyzRead(fin)
    if natoms < params['eigen']:
        params['eigen'] = natoms
    mat = eigenCoulomb(fin,params['eigen'])
    dataxyz.append(mat)
    ids.append(num_files)
    num_files+=1

end = time.time()
tot=end - start
print('para ler '+str(tot))


if opt in range(2,4):
    pdata2 = pd.read_csv(sys.argv[4], sep="\t");
    columns = pdata2.columns.tolist()
    # print(columns)
    pdata2.sort_values(columns[1], axis=0, ascending=True, inplace=True)
    print(pdata2)
    pdata2.sort_values(columns[0], axis=0, ascending=True, inplace=True)
    # ids = pdata2.values[:,0]
    if opt == 2:
        data2 = pdata2.loc[:,['Energy_T']]#.values[:,1:]
    else:
        data2 = pdata2.loc[:,['Energy_T', 'Delta_E ', 'm_T', 'ECN', 'd_av']]#.values[:,1:]
    data2 = data2.values
    X = np.concatenate((data2,np.array(dataxyz)), axis=1)
else:
    X = np.array(dataxyz)

# print(np.shape(X))

# exit()

from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn import metrics
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, silhouette_score


X = StandardScaler().fit_transform(X) 

tsnee = TSNE(n_components=2) ####


X = tsnee.fit_transform(X) ####

kmeans = cluster.KMeans(init='random', #'k-means++', 
        n_clusters=params['n_clusters'], n_init=10, random_state=0)
# kmeans = cluster.KMeans(init='random', 
#         n_clusters=params['n_clusters'])


#########

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
ax1.set_ylim([0, len(X) + (params['n_clusters'] + 1) * 10])



#clusterer =cluster.KMeans(n_clusters=params['n_clusters'], random_state=10) tem o meu
cluster_labels = kmeans.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
silhouette_avg = silhouette_score(X, cluster_labels)
print("For n_clusters =", params['n_clusters'],
          "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
sample_silhouette_values = silhouette_samples(X, cluster_labels)


y_lower = 10
for i in range(params['n_clusters']):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
	ith_cluster_silhouette_values = \
	sample_silhouette_values[cluster_labels == i]

	ith_cluster_silhouette_values.sort()

	size_cluster_i = ith_cluster_silhouette_values.shape[0]
	y_upper = y_lower + size_cluster_i

	color = cm.nipy_spectral(float(i) / params['n_clusters'])
	ax1.fill_betweenx(np.arange(y_lower, y_upper),0, ith_cluster_silhouette_values, facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
	ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
	y_lower = y_upper + 10  # 10 for the 0 samples

ax1.set_title("The silhouette plot for the various clusters.")
ax1.set_xlabel("The silhouette coefficient values")
ax1.set_ylabel("Cluster label")

# The vertical line for average silhouette score of all the values
ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

ax1.set_yticks([])  # Clear the yaxis labels / ticks
ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
colors = cm.nipy_spectral(cluster_labels.astype(float) / params['n_clusters'])
ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7, c=colors, edgecolor='k')

    # Labeling the clusters
centers = kmeans.cluster_centers_
    # Draw white circles at cluster centers
ax2.scatter(centers[:, 0], centers[:, 1], marker='o',c="white", alpha=1, s=200, edgecolor='k')

for i, c in enumerate(centers):
	ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1, s=50, edgecolor='k')

ax2.set_title("The visualization of the clustered data.")
ax2.set_xlabel("Feature space for the 1st feature")
ax2.set_ylabel("Feature space for the 2nd feature")

plt.suptitle(("Silhouette analysis for k-means clustering on sample data of "+str(sys.argv[2])+" " " with n_clusters = %d" % params['n_clusters']), fontsize=14, fontweight='bold')

#plt.show()
plt.savefig('silAl55_'+str(sys.argv[2])+'.png', format='png')
print("Al55")









##################
# Obtaining the most representative molecules
selected1 = []
clusters = []
kmeans.fit(X)
centroids = kmeans.cluster_centers_

if opt==1:
    # for clus in centroids:
    for clus in range(params['n_clusters']):
        idsIn = np.where(kmeans.labels_==clus)
        sel = idsIn[0][0]
        dMin = np.linalg.norm(centroids[clus,:] - X[idsIn[0][0],:])
        for sample in idsIn[0][1:]:
            dist = np.linalg.norm(centroids[clus,:] - X[sample,:])
            if dist < dMin:
                dMin = dist
                sel = sample
        selected1.append(int(ids[sel]))
        clusters.append(idsIn)
else:
    for clus in range(params['n_clusters']):
        idsIn = np.where(kmeans.labels_==clus)
        sel = idsIn[0][0]
        enMin = data2[sel,0]
        for sample in idsIn[0][1:]:
            energy = data2[sample,0]
            if energy < enMin:
                enMin = energy
                sel = sample
        selected1.append(int(ids[sel]))
        clusters.append(idsIn)

print("Selected [op1]: ")
print(selected1)

a=selected1

os.system('mkdir '+str(sys.argv[2])+'/selected_'+str(sys.argv[2]))

for y in a:
    
    print(y)
    yy=baseFolder[y-1]

    print(yy)

    os.system('cp  '+str(yy)+' '+str(sys.argv[2])+'/selected_'+str(sys.argv[2])) #elements may change according to filenames


#end = time.time()
#tot=end - start
#print('tempo total '+str(tot))


exit()
"""
#exit()

# plotting PCA illustration
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
print("oi 5")


#pca = PCA(n_components=2)########################
pca = TSNE(n_components=2)

reduced_X = pca.fit_transform(X)
#print(pca.explained_variance_ratio_)####################### 
# exit() 

kmeans.fit(reduced_X)

print("oi 6")
 
# Step size of the mesh.
h = 0.01
x_min, x_max = reduced_X[:, 0].min() - 1, reduced_X[:, 0].max() + 1
y_min, y_max = reduced_X[:, 1].min() - 1, reduced_X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

print("oi 7")

# Obtain labels for each point in mesh. Use last trained model.

#Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
#Z = Z.reshape(xx.shape)

plt.figure(1)
plt.clf()
print("oi 8")

#plt.imshow(Z, interpolation='nearest',
#           extent=(xx.min(), xx.max(), yy.min(), yy.max()),
#           cmap=plt.cm.Paired,
#           aspect='auto', origin='lower')

plt.plot(reduced_X[:, 0], reduced_X[:, 1], 'k.', markersize=3)
# Plot the centroids as a white X
centroids = kmeans.cluster_centers_

plt.scatter(centroids[:, 0], centroids[:, 1],
            marker='x', s=75, linewidths=3,
            color='b', zorder=10)
print("oi 9")

selected1 = np.array(selected1)-1
for i in selected1:
    plt.text(reduced_X[i, 0], reduced_X[i, 1], str(i+1),fontsize=10)

#y = [2.56422, 3.77284, 3.52623, 3.51468, 3.02199]
#z = [0.15, 0.3, 0.45, 0.6, 0.75]
#n = [58, 651, 393, 203, 123]

#fig, ax = plt.subplots()
#ax.scatter(z, y)

#for i, txt in enumerate(n):
#    ax.annotate(txt, (z[i], y[i]))


plt.title('Clustering Result (PCA-reduced data)\n'
          'Centroids are marked with white cross')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(())
plt.yticks(())

end = time.time()
tot=end - start
print(str(tot))
#plt.show()
plt.savefig('clust'+str(sys.argv[2])+'.pdf', format='pdf')
"""
