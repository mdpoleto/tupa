def grep_efield_value(col, path):
    f = open(path + "/ElecField.dat", 'r')
    lines = f.readlines()
    f.close()

    value = None
    
    for line in lines:
        cols = line.split()
        if cols[0] == "#AVG:":
            value = float(cols[col])
            break
    
    return value


