from cmath import inf
import math
import random
import numpy as np

def G2D(G):
    #相邻栅格点初始化
    #该函数将所有点周围邻域周围路径所需要代价值写入，类似于邻接图标注距离的方式，此处是将所有的格栅的对邻域邻接图都表示出来，图内第i行第j列对于邻域的邻接关系为D[（i-1）*l+j]的一个邻接矩阵来表示
    l=len(G)
    #D是所有点到所有点的矩阵，纵向将所有点作为起点，横向将所有点作为邻域局部到达点
    D=np.zeros((l*l,l*l))
    for i in range(0,l):
        for j in range(0,l):
            if G[i][j] == 0:  #起始栅格是否有障碍
                for m in range(0,l):
                    for n in range(0,l):
                        if G[m][n] == 0: #到达栅格是否有障碍
                            im=abs(i-m)
                            jn=abs(j-n)
                            if(im+jn==1 or (im==1 and jn==1)):
                                D[i*l+j][m*l+n] = (im+jn)**0.5 
                                #将所有点到所有邻域点的代价值填入D中
                                #D[i*l+j][m*l+n]在矩阵中定位邻域矩阵的具体位置
                                #(im+jn)**0.5在具体位置中填入具体代价值 
    return D

f  = open("G.txt")
G = []

for lines in f:
    temp = lines.strip('\n').split(" ")
    for i in range(0,len(temp)) :
        temp[i] = int(temp[i])
    G.append(temp)
f.close()

MM= len(G)
Tau = np.ones((MM*MM,MM*MM))
Tau = 8*Tau
K = 100
M=50                 #蚂蚁个数
S=1                   #最短路径的起始点

E=MM*MM               #最短路径的目的点

Alpha=1               # Alpha 表征信息素重要程度的参数
Beta=7                # Beta 表征启发式因子重要程度的参数
Rho=0.1               #Rho 信息素蒸发系数
Q=1.5                  # Q 信息素增加强度系数 
minkl=inf 
mink=0 
minl=0
D = G2D(G)
N = len(D)           #N表示问题的规模
a = 1                #小方格像素的边长
#计算终止点的横坐标
Ex = a*((E % MM)-0.5)
if Ex ==-0.5:
    Ex = MM-0.5
#计算终止点的纵坐标
Ey = a*(MM+0.5-math.ceil(E/MM))   #其中ceil为向上取整函数
Eta = np.zeros(N)  #启发式信息，取为至目标点的直线距离的倒数
#以下启发式信息矩阵
for i in range(1,N+1):
    ix = a*((i%MM)-0.5)
    if ix == -0.5:
        ix = MM-0.5
    iy = a*(MM+0.5-math.ceil(i/MM))
    if i == E:
        Eta[i-1] = 100
    else:
        Eta[i-1] = 1/math.sqrt((ix-Ex)**2+(iy-Ey)**2)
def print_list(l):
    for i in range(0,20):
        for j in range(0,20):
            print('{:.5f}'.format(l[(i)*20+j]),end=' ')
        print()
def find(l):
    #该函数用于寻找某个列表中的非零值所在列表之中的位置，返回值是一个存储非零值位置的数据列表
    temp = []
    for i in range(0,len(l)):
        if l[i] != 0:
            t = i
            temp.append(t)
    return temp
def find_selecet(l):
    #该函数用于生成一个随机数位于（0,1）之间，用于概率函数的分布选择
    temp=random.random()
    for i in range(0,len(l)):
        if l[i] >temp:
            t = i
            break
    return t
ROUTES =[] #用列表储存每一代的每一只蚂蚁的爬行路线
for i in range(0,K):
    TEMP=[]
    for j in range(0,M):
        TEMP.append([])
    ROUTES.append(TEMP)
PL = np.zeros((K,M))
for k in range(0,K):
    for m in range(0,M):
        #状态初始化
        W = S                 #当前节点初始化为起始点
        Path = []
        Path.append([1,1])              #爬行路线初始化
        PLkm = 0              #爬行路线长度初始化
        TABUkm = np.ones(N)   #禁忌表初始化
        TABUkm[S-1] = 0         #已经在初始点，排除
        DD = D.copy()                #邻接矩阵初始化
        #这里需要利用深拷贝，不能直接赋值，否则下面的内容工会影响该D的本身取值导致程序不能够正常运行
        #下一步可以前往的节点
        DW = DD[W-1]
        DW1 = find(DW)
        for j in DW1:
            #该循环是确认此处的点是否已经存于禁忌表之中，若节点已在禁忌表，则下一个前往节点则不会考虑
            if TABUkm[j] == 0:
                DW[j]=0
        LJD = find(DW)
        Len_LJD = len(LJD)    #可选将诶点的个数
        while ( W+1 !=E and Len_LJD>0 ):
            #转轮赌法选择下一步怎么走
            #PP = np.zeros((Len_LJD))
            PP=[]
            for i in range(0,Len_LJD):
                PT=(Tau[W-1][LJD[i]]**(Alpha))*((Eta[LJD[i]]**(Beta)))
                PP.append(PT)
            sumpp = sum(PP)
            PP = PP/sumpp  #建立概率分布函数
            Pcum=[]
            #Pcum是一个为了便于选择的列表，此种选择不会重复，且最后总概率为1
            Pcum.append(PP[0])
            for i in range(1,Len_LJD):
                Pcum.append(Pcum[i-1]+PP[i])
            Select = find_selecet(Pcum)
            to_visit = LJD[Select]
            Path.append([W+1,to_visit+1])   #增加路径显示
            PLkm = PLkm+DD[W][to_visit] #路径长度增加
            W = to_visit                #蚂蚁移到下一节点
            for kk in range(0,N):
                if TABUkm[kk] == 0:
                    DD[W][kk] = 0
                    DD[kk][W] = 0
            TABUkm[W] = 0  #将已经访问过的节点从禁忌表中删除
            DW=DD[W]
            DW1 = find(DW)
            for j in DW1:
                if TABUkm[j] == 0:
                    DW[j]=0
            LJD = find(DW)
            Len_LJD = len(LJD)    #可选节点的个数
            
            
        
        ROUTES[k][m] = Path
        if Path[-1][-1] == E:
            PL[k][m] = PLkm
            if PLkm <minkl:
                mink = k
                minl = m
                minkl = PLkm
        else:
            PL[k][m] = 0
    #更新信息素：
    Delta_Tau = np.zeros((N,N))  #更新量初始化
    for m in range(0,m):
        if PL[k][m] != 0:
            ROUT = ROUTES[k][m]
            TS = len(ROUT)-1   #跳数（还不知道啥意思）
            PL_km = PL[k][m]
            for s in range(0,TS):
                x = ROUT[s][0]
                y = ROUT[s][1]
                Delta_Tau[x][y] = Delta_Tau[x][y]+ Q/PL_km
                #Delta_Tau[y][x] = Delta_Tau[y][x]+ Q/PL_km
    Tau = (1-Rho)*Tau+Delta_Tau #挥发一部分，增加一部分



print("mink={}".format(mink))
print("minl={}".format(minl))
print("minkl={}".format(minkl))
print(ROUTES[mink][minl])
G_2 = G.copy()
file = open("zhexian.txt",'w')
for pr in ROUTES[mink][minl]:
    px = int((pr[1]-0.5)%20+0.5)
    py = int(pr[1]/20)+1
    if py != 21:
        file.write("{},{}".format(px,py))
        file.write("\n")
file.close()
file = open("lujing.txt",'w')
for i in range(0,20):
    for j in range(0,20):
        if j != 19:
            str_w = str(G_2[i][j])+' '
        else:
            str_w = str(G_2[i][j])
        file.write("{}".format(str_w))
    file.write("\n")
file.close()

