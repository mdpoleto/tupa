import numpy as np


def add_to_dict(dict, key, array):
    if key in dict.keys():
        prev_array = dict[key]
        new_array = np.vstack((prev_array, array))
        dict[key] = new_array
    else:
        dict[key] = array

    return dict

def mag(vector):
    """ Returns the magnitude of the vector.  """
    mag = np.linalg.norm(vector)
    return mag


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def alignment(value1,value2):
    """ Returns % alignment between 2 values """
    # V2 is the total reference
    aligned = abs(value1/value2)*100
    return aligned


def projection(v1,v2):
    """ Returns projection of vector v1 onto vector v2 """
    #projection_u_on_v = (np.dot(u, v)/np.dot(v, v))*v
    proj = v2*(np.dot(v1, v2)/np.dot(v2,v2))
    return proj


array = np.array([
[ 2.7928,  4.5620  , -3.6690],
[-3.3893, -5.5731  ,  4.4337],
[ 0.5705,  0.9162  , -0.8287],
[ 0.5017,  0.9196  , -0.7135],
[ 0.6051,  0.8867  , -0.7044],
[ 4.3292,  7.2801  , -5.4760],
[-3.6098, -6.0334  ,  4.5394],
[-0.0303, -0.0603  ,  0.0447],
[ 4.3498,  6.2199  , -5.4316],
[-3.1755, -4.5356  ,  3.9678],
[ 3.0902,  3.8129  , -3.6537],
[-3.0968, -3.8158  ,  3.6632],
[-0.7049, -0.8854  ,  0.8292],
[-0.5485, -0.6636  ,  0.6516],
[ 4.6819,  7.6133  , -4.7239],
[-4.8002, -7.8140  ,  4.8436],
[ 0.0646,  0.1027  , -0.0622],
[ 0.0653,  0.1235  , -0.0695],
[ 5.2218,  7.1779  , -4.2855],
[-5.3377, -7.3538  ,  4.3806],
[ 0.2276,  0.2805  , -0.1950],
[ 0.2154,  0.2857  , -0.1505],
[ 8.1164,  11.8534 , -6.4571],
[-8.0711, -11.8625 ,  6.4948],
[ 5.6676,  9.6703  , -3.8635],
[-6.1299, -10.4317 ,  4.1726],
[ 0.2280,  0.3523  , -0.1296],
[ 0.1895,  0.3689  , -0.1514],
[ 0.2469,  0.4413  , -0.1449],
[-0.3168, -0.5096  ,  0.2866],
[-0.3656, -0.4851  ,  0.2580],  ])

tmp_dict_res = {"MET_1": array}


for r, contribution in tmp_dict_res.items():
    resEf  = np.sum(contribution, axis=0)
    resEfmag = mag(resEf)
    # calculate the projection and alignment for each residue
    #resEfproj = projection(resEf,totalEf) # resEf can have a higher magnitude than the total field.
    #resEfprojmag = mag(resEfproj)
    #resEfalignment = alignment(resEfprojmag,totalEfmag) # resEf can have a higher magnitude than the total field,
                                                          # which means % can be > 100%


    print(resEf)
