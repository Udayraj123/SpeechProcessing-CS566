from mpl_toolkits import mplot3d # for projection='3d' to work
from matplotlib import pyplot as plt
import pandas as pd
from scipy.io import loadmat

# fig = plt.figure()
# plotNum=1
def plotEmbeddings(df):
	# global fig
	# global plotNum
	df = df.drop(df.columns[0],axis=1);
	# print(df)
	x,y,z=df[1],df[2],df[3]
	# ax = fig.add_subplot(1, 2, plotNum, projection='3d')
	# plotNum+=1	
	ax = plt.axes(projection='3d')
	ax.scatter3D(x,y,z,cmap='Greens')
	# ax.scatter(x,y,cmap='Greens')
	plt.show()

colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'];
markers=[".",",","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h","H","+","x","X","D","d","|","_"];

def plotEmbeddingsWithColors(df,labelDict,title="3D embeddings"):	
	df_=pd.concat([df,labelDict],axis=1)
	last=df_.columns.size-1	
	df_.columns=range(last+1)
	# print(df_)
	# markers=['o','^','o','^','o','^','o'];
	# print(df_[last])
	# df_[last-1]=df_[last-1].apply(lambda x:markers[x])
	df_[last]=df_[last].apply(lambda x:colors[x])
	# nodes 

	x,y,z,c=df[1],df[2],df[3],df_[last]	

	ax = plt.axes(projection='3d')
	ax.scatter3D(x,y,z,c=c)
	ax.set_title(title)
	plt.savefig(title+".png")
	plt.show()

def plotEmbeddingsMultiClass(df,labelDict,title="3D embeddings"):	
	start=df.columns.size	
	df_=pd.concat([df,labelDict],axis=1)
	last=df_.columns.size-1	
	df_.columns=range(last+1)
	colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'];
	def getColor(row):
		d=1
		num=0
		# print(row[start:last+1])
		for i in range(start,last+1):
			num += row[i]*d
			d*=2
		num=int(num)
		return colors[num % len(colors)]

	df_[last]=df_.apply(getColor, axis=1)

	x,y,z,c=df[1],df[2],df[3],df_[last]	

	ax = plt.axes(projection='3d')
	ax.scatter3D(x,y,z)#,c=c)
	ax.set_title(title)

	plt.savefig(title+".png")
	plt.show()

def setIndex(df):
	# Assuming node ids go from 0 to N-1
	df = df.set_index(df.columns[0])
	df.sort_index(inplace=True)
	return df

# df1 = pd.read_csv("deepwalk/emb/homo_sapiens3dmat_dw.emb",delimiter=' ',skiprows=1,header=None)
df1 = pd.read_csv("deepwalk/emb/homo_sapiens3d_dw.emb",delimiter=' ',skiprows=1,header=None)
df1 = setIndex(df1)
df2 = pd.read_csv("node2vec/emb/homo_sapiens3d_n2v.emb",delimiter=' ',skiprows=1,header=None)
df2 = setIndex(df2)
labelDict1 = pd.read_csv("data/homo_sapiens_labels.txt",delimiter=' ',header=None)
labelDict1=labelDict1.drop(labelDict1.columns[0],axis=1)
plotEmbeddingsMultiClass(df1,labelDict1,title="3D embeddings for PPI dataset using Node2Vec")
plotEmbeddingsMultiClass(df2,labelDict1,title="3D embeddings for PPI dataset using Deepwalk")

df1 = pd.read_csv("deepwalk/emb/cora3d_dw.emb",delimiter=' ',skiprows=1,header=None)
df1 = setIndex(df1)
df2 = pd.read_csv("node2vec/emb/cora3d_n2v.emb",delimiter=' ',skiprows=1,header=None)
df2 = setIndex(df2)
labelDict2 = pd.read_csv("data/cora_labels.txt",delimiter=' ',header=None)
labelDict2=labelDict2[[1]]
plotEmbeddingsWithColors(df1,labelDict2,title="3D embeddings for Cora dataset using Node2Vec")
plotEmbeddingsWithColors(df2,labelDict2,title="3D embeddings for Cora dataset using Deepwalk")
