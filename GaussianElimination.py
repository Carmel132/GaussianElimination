'''
3 elementary row operations

swap 2 rows
multiply a row by a nonzero factor
add the sum of a row onto another row
'''
from enum import Enum
from copy import deepcopy
class Operations(Enum):
    SWAP=0
    SCALE=1
    ADD=2

class Operation:
    def __init__(_, type, x1, x2, x3=None):
        _.type = type
        _.x1 = x1
        _.x2 = x2
        _.x3 = x3
    def __repr__(_):
        return f"[{_.type}, {_.x1}, {_.x2}, {_.x3}]"
class Row:
    def __init__(_, _vals:list[float]):
        _.vals = _vals
    def __mul__(_, k:float):
        return Row([x * k for x in _.vals])
    def __truediv__(_, k:float):
        return Row([x / k for x in _.vals])
    def __add__(_, right):
        return Row([i + j for i, j in zip(_.vals, right.vals)])
    def __sub__(_, right):
        return Row([i - j for i, j in zip(_.vals, right.vals)])
    def __repr__(_):
        return f"[{', '.join(map(str, _.vals))}]"
    def __str__(_):
        return f"[{', '.join(map(str, _.vals))}]"
    def __getitem__(_, i):
        return _.vals[i]
    def __setitem__(_, i, v):
        _.vals[i] = v
    def __eq__(_, other):
        return all(map(lambda x: x[0] == x[1], zip(_.vals, other.vals)))
    def __ne__(_, other):
        return any(map(lambda x: x[0] != x[1]), zip(_.vals, other))
    def pivotIdx(_):
        for i, k in enumerate(_.vals):
            if k != 0: return i
        return -1
    def pivot(_):
        for i in _.vals:
            if i != 0:return i
        return 0
    def isZero(_):
        return not any(_.vals)
    
    

class Matrix:
    def __init__(_, height:int, length:int):
        _.height = height
        _.length = length
        _.vals = [Row([0] * length) for i in range(height)]
        _.operations = []
    def __repr__(_):
        ret = []
        for i, row in enumerate(_.vals):
            r = []
            for j, v in enumerate(row.vals):
                maxCharacters = len(str(max(_.getColumn(j).vals, key=lambda x:len(str(x)))))
                r.append(" " * (maxCharacters - len(str(v))) + str(v))
            ret.append( f"[{', '.join(r)}]" )
        return "\n".join(ret)
    def swap(_, i1, i2):
        _.vals[i1], _.vals[i2] = _.vals[i2], _.vals[i1]
        _.operations.append(Operation(Operations.SWAP, i1, i2))
    def scale(_, i, k):
        _.vals[i] = _.vals[i] * k
        _.operations.append(Operation(Operations.SCALE, i, k))
    def add(_, i1, i2, f=1): #adds i2 * f onto i1
        _.vals[i1] = _.vals[i1] + (_.vals[i2] * f)
        _.operations.append(Operation(Operations.ADD, i1, i2, f))
    def setRow(_, i, row):
        _.vals[i] = row
    def getRow(_, i):
        return _.vals[i]
    def order(_):
        _.vals = sorted(_.vals, key=lambda x: _.height + 1 if x.isZero() else x.pivotIdx())
    def RREF(_):
        _.REF()
        for i in range(min(_.height, _.length)):
            r = _.getRow(i)
            p = r.pivotIdx()
            _.setRow(i, r / r.pivot())
            for j in range(_.height):
                if i == j: continue
                _.add(j, i, -_.getRow(j)[p])
    def REF(_):
        h = 0
        k = 0
        while h < _.height and k < _.length:
            i_max = max(range(h, _.height), key=lambda x: abs(_.vals[x][k]))
            if _.vals[i_max][k] == 0: k = k+1
            else:
                _.swap(h, i_max)
                for i in range(h + 1, _.height):
                    f = _.vals[i][k] / _.vals[h][k]
                    _.vals[i][k] = 0
                    for j in range(k + 1, _.length):
                        _.vals[i][j] = _.vals[i][j] - _.vals[h][j] * f
                h = h + 1
                k = k + 1
    @staticmethod 
    def identity(n):
        ret = Matrix(n, n)
        for i in range(n):
            ret.setRow(i, Row([1.0 if i == j else 0.0 for j in range(n)]))
        return ret
    def operation(_, o:Operation):
        if o.type == Operations.SWAP:
            _.swap(o.x1, o.x2)
        if o.type == Operations.SCALE:
            _.scale(o.x1, o.x2)
        if o.type == Operations.ADD:
            _.add(o.x1, o.x2, o.x3)
    def getColumn(_, i):
        return Row([r[i] for r in _.vals])
    def getColumns(_):
        return [_.getColumn(i) for i in range(_.length)]
    def multiply(_, other):
        ret = Matrix(_.height, other.length)
        for i, r in enumerate(_.vals):
            for j, c in enumerate(other.getColumns()):
                ret.vals[i].vals[j] = dot(r, c)
        return ret
    def getE(_, iterations=5):
        c = deepcopy(_)
        for i in range(iterations):
            c.RREF() #hopefully this works lmfao
        ret = Matrix.identity(c.height)
        for o in reversed(c.operations):
            other = Matrix.identity(c.height)
            other.operation(o)
            ret = ret.multiply(other)
        return ret
    def getInverse(_, iterations=5):
        if not _.isInversible(iterations):return -1
        c = deepcopy(_)
        for i in range(iterations):
            c.RREF()
        ret = Matrix.identity(c.height)
        for o in c.operations:
            ret.operation(o)
        return ret
    def isInversible(_, iterations = 5):
        if _.height != _.length: return False
        c = deepcopy(_)
        for i in range(iterations):
            c.RREF()
        return c.vals == Matrix.identity(c.height).vals
    def round(_): # rounds to nearest 10^-16
        for r in range(_.height):
            for v in range(_.length):
                _.vals[r][v] = round(_.vals[r][v], 16)
    def getTransverse(_):
        ret = Matrix(_.length, _.height)
        for i in range(_.length):
            c = _.getColumn(i)
            c.vals = list(reversed(c.vals))
            ret.setRow(i, c)
        return ret
    def getRightInverse(_):
        return _.getTransverse().multiply((_.multiply(_.getTransverse())).getInverse())
    def getLeftInverse(_):
        return (_.getTransverse().multiply(_)).getInverse().multiply(_.getTransverse())
def dot(v1 ,v2): return sum(map(lambda x: x[0] * x[1], zip(v1, v2)))

a = Matrix(3, 2)
a.setRow(0, Row([1, 2]))
a.setRow(1, Row([0, -1]))
a.setRow(2, Row([1, 1]))
a.RREF()
#al = a.getLeftInverse()
print(a)
#print("\n", al.multiply(a))