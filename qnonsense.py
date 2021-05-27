import math
from sympy.solvers import solve
from sympy import Symbol

#basically unrunnable if you don't have a beefier computer than i do 
#but i'm putting it in the repo anyway
def pmaker(data, alpha, valuelist):
    """intended to do our fake binary p(x) builder. does it work? who knows"""
    qlist = []
    x = Symbol('x')
    norm_scriptx = len(valuelist)**(len(data))
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))
    for value in valuelist:
        mg=0
        for i in range(0, len(data)-data.count(value)+1):
            mg += math.comb(len(data),i)*((len(valuelist)-1)**i)
        mg*len(valuelist)
        nu = norm_scriptx/mg
        s_lowerbound = alpha*nu/(r*u)
        q = solve(x**data.count(value)*(1-x)**(len(data)-data.count(value)) - s_lowerbound*u, x)
        qlist += [q[0]]
    return qlist



