from sage.calculus.integration import numerical_integral

def dimMapping(order):
    mapping = { 1                      : str(1), 
                2                      : str(2), 
                order                   : str(order) ,
                order+1                 : str(order+1) ,
                (order+1)*(order+2)/2   : str((order+1)*(order+2)/2),
                (order+1)*(order+2)/2+1 : str((order+1)*(order+2)/2+1),
                (order+1)*(order+2)     : str((order+1)*(order+2)),
                (2*order+1)             : str((2*order+1)),
                (2*order+1)**2          : str((2*order+1)**2)}
    return mapping[order]
    
def numerical_integral_2d(func,a1=0,b1=1,a2=0,b2=1,c=-1):
    integrand = lambda x : numerical_integral(func,a2,b2+x*c,params=[x])[0]
    return numerical_integral(integrand,a1,b1)[0]

def generatePatch(order,show):
    num = 2*order+1
    nodes = []
    for i in range(0,num+1):
        nodes.append([])
        for j in range(0,i+1):
            nodes[i].append([(i-j)/float(num),j/float(num)])
    triangles = [[[0,0],[0,0],[0,0]]] * num * num
    
    ind = 0
    for i in range(0,num):
        nodes_top    = nodes[i]
        nodes_bottom = nodes[i+1]
        for j in range(0,i+1):
            if not j == 0:
                ind = ind + 2
            triangles[ind] = [nodes_top[j], nodes_bottom[j], nodes_bottom[j+1]]
        ind = ind + 1    
            
    ind = 2
    for i in range(1,num):
        nodes_top    = nodes[i]
        nodes_bottom = nodes[i+1]
        for j in range(1,i+1):
            if not j == 1:
                ind = ind + 2
            triangles[ind] = [nodes_top[j], nodes_top[j-1], nodes_bottom[j]]
        ind = ind + 3         
    i = 0
    if show:
        for triangle in triangles:
            plt.plot([triangle[0][0], triangle[1][0], triangle[2][0], triangle[0][0]],
                    [triangle[0][1], triangle[1][1], triangle[2][1], triangle[0][1]])
            center_x = (triangle[0][0]+ triangle[1][0] + triangle[2][0])/3.0
            center_y = (triangle[0][1]+ triangle[1][1] + triangle[2][1])/3.0
            plt.text(center_x,center_y,str(i))
            i = i + 1 

        plt.show()
    return triangles

def getFVEdgeIndeces(order,edge):
    indeces = []
    num = 2 * order + 2 
    if edge == "l":
        for i in range(2*order,-1,-1):
            indeces.append(i*i)
    if edge == "m":
        for i in range(2*order,-1,-1):
            indeces.append((2 * order + 1) **2 - (2 * order + 1) - (2 * order ) + 1 * i * 2)
    if edge == "r":
        for i in range(0,2*order+1):
            indeces.append((i+1)*(i+1)-1)
    return indeces
getFVEdgeIndeces(3,"r")

def getDGIndeces(order,direction):
    indeces=[]
    if direction == "l":
        return list(range(0,order + 1))
    if direction == "m":
        idx=order
        for i in range(0,order + 1):
            indeces.append(idx)
            idx = idx + order - i
    if direction == "r":
        idx=0
        for i in range(0,order + 1):
            indeces.append(idx)
            idx = idx + order + 1 - i
    return indeces
        
