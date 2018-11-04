from mpl_toolkits import mplot3d # for projection='3d' to work
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = (15.0, 10.0)
import time
import pandas as pd
import numpy as np
from matplotlib.lines import Line2D

n_components=2
save=False
vowels=['a', 'e', 'i', 'o', 'u'];
colors=['r', 'g', 'b', 'y', 'm'];
legend_elements = [Line2D([0], [0], marker='o', color='w', label=v, markerfacecolor=c, markersize=15) for v,c in zip(vowels,colors)]
from sklearn import decomposition
from sklearn.manifold import TSNE
def plotTSNE(ax,df,labelCols, init = 'random',perplexity=50, n_components=2, n_iter=2000):	
	if(df.shape[0]!=labelCols.shape[0]):
		print("Incompatible shapes!", df.shape, labelCols.shape)
		return
	# concat the label col(s)
	df_=pd.concat([df,labelCols],axis=1)
	last=df_.columns.size-1
	df_.columns=range(last+1)
	df_[last]=df_[last].apply(lambda x:colors[x])
	
	
	t_prev = time.time()
	# learning_rate n_iter=1000
	print("tSNE running...","Perplexity :",perplexity, "Shape:")
	tsne_obj = TSNE(n_components=n_components,n_iter=n_iter, random_state=0, init = init, perplexity=perplexity)
	print(df.shape)
	tsneVecs = tsne_obj.fit_transform(df)
	tsneVecs = tsneVecs.T 
	print("Iterations: ",tsne_obj.n_iter_, "KL Divergence: %.2f" % (tsne_obj.kl_divergence_) )
	print("Time taken: ", "%.2f"%(time.time()-t_prev))
	print("Plotting subplot...")	
	if(n_components==2):
		# ax.scatter(X[0],X[1],c=df_[last])
		ax.scatter(tsneVecs[0],tsneVecs[1],c=df_[last])
	    # ax.set_xlabel("X")
	    # ax.set_ylabel("Y")	
	else:
		ax = plt.axes(projection='3d')
		ax.scatter3D(tsneVecs[0],tsneVecs[1],tsneVecs[2],c=df_[last])
	ax.legend(handles=legend_elements)
	ax.set_title("tSNE Plot with Perplexity: "+str(perplexity))

df2 = pd.read_csv("coeffs/universe.txt",delimiter='\t',header=None)
df2 = df2.drop(df2.columns[0],axis=1);
df2 = df2.dropna(axis=1)
labels = pd.DataFrame(np.array([[i]*50 for i in range(5)]).flatten())
# print(labels)
fig, ax = plt.subplots(2, 2)
for i,p in enumerate([15,30,40,50]):
	plotTSNE(ax[i//2,i%2],df2,labels,n_components=n_components, init = 'pca', perplexity=p)
if(save):
	print("Saving..")
	plt.savefig("tsne_vowel_plots_lowp.png")
else:
	plt.show()

fig, ax = plt.subplots(2, 2)
for i,p in enumerate([80,100,120,160]):
	plotTSNE(ax[i//2,i%2],df2,labels,n_components=n_components,init = 'pca', perplexity=p)

if(save):
	print("Saving..")
	plt.savefig("tsne_vowel_plots_highp.png")
else:
	plt.show()
