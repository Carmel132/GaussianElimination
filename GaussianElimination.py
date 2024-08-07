'''
3 elementary row operations

swap 2 rows
multiply a row by a nonzero factor
add the sum of a row onto another row
'''
from enum import Enum
from copy import deepcopy
from math import sqrt, floor, pow, log10
from functools import reduce
class Frac:
    def __init__(_, n, d=1, neg = False):
        _.nu = int(abs(n))
        _.de = int(abs(d))
        _.negative = (n < 0 or d < 0) or neg
        _.isZero = lambda: _.nu == 0
        _.isInfinity = lambda: _.de == 0
    def smpy(_):
        if _.isZero(): _.de = 1
        else:
            i = min(_.nu, _.de)
            while 2 <= i:
                if _.nu % i == 0 and _.de % i == 0:
                    _.nu = int(_.nu / i)
                    _.de = int(_.de / i)
                i -= 1
        return _
    @staticmethod
    def fromFloat(f):
        if type(f) is Frac: return f
        count = 0
        while f - floor(f) != 0:
            f *= 10
            count += 1
        return Frac(f, pow(10, count)).smpy()
    def flip(_):
        return Frac(_.de, _.nu, _.negative)
    def __repr__(_):
        return f"{'-' if _.negative and not _.isZero() else ''}{_.nu}{f'/{_.de}' if _.de != 1 else ''}"
    def __str__(_):
        return f"{'-' if _.negative and not _.isZero() else ''}{_.nu}{f'/{_.de}' if _.de != 1 else ''}"
    def __mul__(_, o):
        return Frac(_.nu * o.nu, _.de * o.de, _.negative ^ o.negative).smpy()
    def __truediv__(_, o):
        return Frac(_.nu * o.de, _.de * o.nu, _.negative ^ o.negative).smpy()
    def __add__(_, o):
        return Frac(_.nu * o.de * (-1 if _.negative else 1) + o.nu * _.de * (-1 if o.negative else 1), _.de * o.de).smpy()
    def __sub__(_, o):
        return Frac(_.nu * o.de * (-1 if _.negative else 1) - o.nu * _.de * (-1 if o.negative else 1), _.de * o.de).smpy()
    def __eq__(_, o):
        return (_ - o).isZero()
    def __neq__(_, o):
        return not (_ - o).isZero()
    def __lt__(_, o):
        return (_ - o).negative and _ != o
    def __gt__(_, o):
        return not (_ - o).negative and _ != o
    def __bool__(_):
        return not _.isZero()
    def __abs__(_):
        return Frac(_.nu, _.de)
    def __float__(_):
        return _.nu / _.de * (-1 if _.negative else 1)
    def __neg__(_):
        return Frac(_.nu, _.de, not _.negative)

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
        _.vals = list(map(Frac.fromFloat, _vals))
    def __mul__(_, k:Frac):
        return Row([x * k for x in _.vals])
    def __truediv__(_, k:Frac):
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
    def __len__(_):
        return len(_.vals)
    def pivotIdx(_):
        for i, k in enumerate(_.vals):
            if not k.isZero(): return i
        return -1
    def pivot(_):
        for i in _.vals:
            if not i.isZero(): return i
        return Frac(0)
    def isZero(_):
        return not any(_.vals)
    def __neg__(_):
        return Row(map(lambda x: -x, _.vals))

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
    @staticmethod
    def fromRows(rs):
        ret = Matrix(len(rs), len(rs[0]))
        for i, r in enumerate(rs):
            ret.setRow(i, r)
        return ret
    def order(_):
        _.vals = sorted(_.vals, key=lambda x: _.height + 1 if x.isZero() else x.pivotIdx())
    def RREF(_):
        _.EF()
        for i in range(min(_.height, _.length)):
            r = _.getRow(i)
            if r.isZero():continue
            p = r.pivotIdx()
            _.scale(i, r.pivot().flip())
            for j in range(_.height):
                if i == j: continue
                _.add(j, i, -_.getRow(j)[p])
    def EF(_):
        h = 0
        k = 0
        while h < _.height and k < _.length:
            i_max = max(range(h, _.height), key=lambda x: abs((_.vals[x][k])))
            if _.vals[i_max][k].isZero(): 
                k = k+1
            else:
                _.swap(h, i_max)
                for i in range(h + 1, _.height):
                    _.add(i, h,  -_.vals[i][k] / _.vals[h][k])
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
    def getE(_):
        c = deepcopy(_)
        c.RREF() #hopefully this works lmfao
        ret = Matrix.identity(c.height)
        for o in reversed(c.operations):
            other = Matrix.identity(c.height)
            other.operation(o)
            ret = ret.multiply(other)
        return ret
    def E(_):
        ret = Matrix.identity(_.height)
        for o in reversed(_.operations):
            other = Matrix.identity(_.height)
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
    def getTranspose(_):
        ret = Matrix(_.length, _.height)
        for i in range(_.length):
            c = _.getColumn(i)
            c.vals = list(reversed(c.vals))
            ret.setRow(i, c)
        return ret
    def getRightInverse(_):
        return _.getTranspose().multiply((_.multiply(_.getTranspose())).getInverse())
    def getLeftInverse(_):
        return (_.getTranspose().multiply(_)).getInverse().multiply(_.getTranspose())
    def getRowSpaceBasis(_):
        ret = []
        c = deepcopy(_)
        c.RREF()
        for r in c.vals:
            if not r.isZero(): ret.append(r)
        return ret
    def getColumnSpaceBasis(_):
        c = deepcopy(_)
        c = c.getTranspose().getTranspose().getTranspose()
        c.RREF()
        ret = []
        for r in c.vals:
            if not r.isZero(): ret.append(r)
        return ret
    def getNullTransposeSpaceBasis(_):
        c = deepcopy(_)
        c.EF()
        e = c.E()
        ret = []
        for i, r in enumerate(c.vals):
            if r.isZero(): ret.append(e.getRow(i))
        return ret
    def getNullSpaceBasis(_):
        c = deepcopy(_)
        c.RREF()
        def isPivotColumn(c:Row):
            c1 = c.vals.count(Frac(1))
            c0 = c.vals.count(Frac(0))
            return c1 == 1 and c1 + c0 == len(c.vals)
        freeVars = [i for i in range(c.length) if not isPivotColumn(c.getColumn(i))]
        ret = Matrix(c.length, c.length)
        idx = 0
        for i in range(c.length):
            if i in freeVars:
                ret.getRow(i)[i] = 1
            else:
                row = c.getRow(idx)
                row[i] = 0
                ret.setRow(i, -row)
                idx += 1
        return list(filter(lambda x: not x.isZero(), ret.getColumns()))

def dot(v1 ,v2): return sum(map(lambda x: x[0] * x[1], zip(v1, v2)), start=Frac(0))

a = Matrix.fromRows([
    Row([-3, 6, -1, 1, -7]),
    Row([1, -2, 2, 3, -1]),
    Row([2, -4, 5, 8, -4])
    ])
#print(a)
print(a.getNullSpaceBasis())
#a.RREF()
#print(a)
#print("A:", a, sep="\n")
#print("Row:", Matrix.fromRows(a.getRowSpaceBasis()).getTranspose(), sep="\n")
#print("Column:", Matrix.fromRows(a.getColumnSpaceBasis()).getTranspose(), sep="\n")
#print("Null:", Matrix.fromRows(a.getNullSpaceBasis()).getTranspose(), sep="\n")
#print("Null Transpose:", Matrix.fromRows(a.getNullTransposeSpaceBasis()).getTranspose(), sep="\n")