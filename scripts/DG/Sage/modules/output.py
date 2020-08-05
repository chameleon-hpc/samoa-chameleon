
def functionToFrotranString(expr,i):
    basis_name = "basis_"+str(i)
    results_str  ="function "+basis_name+"(x,y)\n"
    results_str +="  real(kind = GRID_SR) ::"+basis_name +"\n"
    results_str +="  real(kind = GRID_SR),intent(in) :: x\n"
    results_str +="  real(kind = GRID_SR),intent(in) :: y\n"
    results_str +="  "+basis_name+" = "

    for coeff_y, exp_y in expr.coefficients(y):
        if exp_y != 0:
            results_str += "y**"+str(exp_y) + " * "
        results_str += "("
        for coeff_x, exp_x in coeff_y.coefficients(x):
            if exp_x != 0:
                results_str += "x**"+str(exp_x)
                results_str += " * "
            results_str += str(coeff_x)+"_GRID_SR + "
        results_str = results_str[:-3] # remove last + for x coeff

        results_str += ") + "
    results_str  = results_str[:-3] # remove last + for y coeff
    results_str +="\nend function "+basis_name

    return results_str

def printBasisFunctions(basis2D,order):
    basis_functions_string = ""
    for basis_idx in range(0,len(basis2D)):
        expr=basis2D[basis_idx]
        basis_functions_string += functionToFrotranString(expr,basis_idx+1)
        basis_functions_string += "\n \n"

    basis_functions_string +=\
"""function evaluate(Q,coords)
  real(kind = GRID_SR)            :: evaluate
  real(kind = GRID_SR),intent(in) :: Q(_SWE_DG_DOFS)
  real(kind = GRID_SR),intent(in) :: coords(2) 

  evaluate =&
  """
    for basis_idx in range(0,len(basis2D)):
        basis_functions_string += "basis_"+str(basis_idx+1)+"(coords(1),coords(2))"
        basis_functions_string += " * Q("+str(basis_idx+1)+") +&"
        basis_functions_string +="\n  "
    
    basis_functions_string = basis_functions_string[:-6]
    basis_functions_string += "\nend function evaluate" 

    with open("basis_"+str(order)+".incl","w") as output:
        output.write(basis_functions_string)

def getTensorShape(tensor,rank):
    shape=()
    tensor_slice=tensor
    for i in range(rank):
        shape= shape[0:i+1] + (len(tensor_slice),)
        tensor_slice=tensor_slice[0]
    return shape

def getJSONFormatTensor(tensor,rank,name):
    epsilon=10.0e-15
    shape = list(getTensorShape(tensor,rank))
    num_entries = 1
    #span index space
    indeces = [(index,) for index in range(0,shape[0])]
    for shape_range in shape[1:]:
        new_indeces=[]
        for index in indeces:
            for new_index in range(0,shape_range):
                new_indeces += [index + (new_index,)]
        indeces=new_indeces
        #print(indeces)
    
    string="{\n"
    string+='  "name": "{}",\n'.format(name)
    string+='  "rank":  {},\n' .format(rank)
    string+='  "shape": {} ,\n'.format(shape)
    string+='  "entries": ['
    for index in indeces:
        entry=tensor[index[0]]
        for i in range(1,len(index)):
            entry=entry[index[i]]
        if(abs(entry)>epsilon):
            string += "["
            for i in range(0,len(index)):
                string += '{}, '.format(index[i]+1) 
            string += '"{}"'.format(entry)
            string += "],"

        #string+= '{}, {}, "{}"'.format(row+1,column+1,entry)
    #        string+=']'
    string=string[:-1]#remove last comma
    string+=']\n'
    string+="}"
    return string

def getJSONFormat(matrix,name):
    epsilon=10.0e-15
    string="{\n"
    string+='  "name": "{}",\n'.format(name)
    string+='  "rows":  {},\n'.format(len(matrix.rows()))
    string+='  "columns": {},\n'.format(len(matrix.columns()))
    string+='  "entries": ['
    for row in range(0,len(matrix.rows())):
        if (abs(matrix[row,:].numpy()) > epsilon).any():
            for column in range(0,len(matrix.columns())):
                entry=matrix[row,column]
                if(abs(entry) > epsilon):
                    string+= '[{}, {}, "{}"],'.format(row+1,column+1,entry)
    string=string[:-1]#remove last comma
    string+=']\n'
    string+="}\n"
    return string

def printVector(vector,name,order):
    dimString = dimMapping[vector.degree()]
    vector_name = name+"("+dimString+")"
    output_string="real(kind=GRID_SR),Parameter :: "+vector_name+" = (/ &\n"
    for i in range(vector.degree()):
            #if i != 0 and i % 10 == 0:
            #        output_string = output_string + "&\n"
        output_string = output_string + str(vector[i]) + "_GRID_SR ,"
    output_string = output_string + "&\n"
    #remove last seperator
    output_string = output_string[:-3]
    output_string = output_string + "  /)"
    with open(str(name)+"_"+str(order)+".incl","w") as output:
        output.write(output_string)


def functionToFrotranString(expr,i):
    basis_name = "basis_"+str(i)
    results_str  ="function "+basis_name+"(x,y)\n"
    results_str +="  real(kind = GRID_SR) ::"+basis_name +"\n"
    results_str +="  real(kind = GRID_SR),intent(in) :: x\n"
    results_str +="  real(kind = GRID_SR),intent(in) :: y\n"
    results_str +="  "+basis_name+" = "

    for coeff_y, exp_y in expr.coefficients(y):
        if exp_y != 0:
            results_str += "y**"+str(exp_y) + " * "
        results_str += "("
        for coeff_x, exp_x in coeff_y.coefficients(x):
            if exp_x != 0:
                results_str += "x**"+str(exp_x)
                results_str += " * "
            results_str += str(coeff_x)+"_GRID_SR + "
        results_str = results_str[:-3] # remove last + for x coeff

        results_str += ") + "
    results_str  = results_str[:-3] # remove last + for y coeff
    results_str +="\nend function "+basis_name

    return results_str

def printBasisFunctions(basis2D,order):
    basis_functions_string = ""
    for basis_idx in range(0,len(basis2D)):
        expr=basis2D[basis_idx]
        basis_functions_string += functionToFrotranString(expr,basis_idx+1)
        basis_functions_string += "\n \n"

    basis_functions_string +=\
"""function evaluate(Q,coords)
  real(kind = GRID_SR)            :: evaluate
  real(kind = GRID_SR),intent(in) :: Q(_SWE_DG_DOFS)
  real(kind = GRID_SR),intent(in) :: coords(2) 

  evaluate =&
  """
    for basis_idx in range(0,len(basis2D)):
        basis_functions_string += "basis_"+str(basis_idx+1)+"(coords(1),coords(2))"
        basis_functions_string += " * Q("+str(basis_idx+1)+") +&"
        basis_functions_string +="\n  "
    
    basis_functions_string = basis_functions_string[:-6]
    basis_functions_string += "\nend function evaluate" 

    with open("basis_"+str(order)+".incl","w") as output:
        output.write(basis_functions_string)
        
