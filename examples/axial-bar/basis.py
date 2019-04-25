def N(xi, i, enum):
    if enum == 1:
        if i == 1:    
            return (2*xi-1)**2
        if i == 2:
            return -2*xi*(3*xi-2)
        if i == 3:
            return 2*xi*xi        

    if enum == 2:
        if i == 1:
            return (2*xi-2)*(xi-1)
        if i == 2:
            return -6*xi**2 + 8*xi - 2
        if i == 3:
            return (2*xi-1)**2

def Nprime(xi, i, enum):
    if enum == 1:
        if i == 1:    
            return 4*(2*xi-1)
        if i == 2:
            return 4*(1-3*xi)
        if i == 3:
            return 4*xi

    if enum == 2:
        if i == 1:
            return 4*(xi-1)
        if i == 2:
            return 4*(2-3*xi)
        if i == 3:
            return 4*(2*xi-1)

def H(xi, L, i):
    if i == 2:
        return xi/L
    if i == 1:
        return 1.0 - xi/L
    else:
        raise

def Hprime(xi, L, i):
    if i == 2:
        return 1.0/L 
    if i == 1:
        return -1.0/L


