class Derivada:
    def __init__(self, f, metodo ="adelante", dx= 0.001):
        self.f = f
        self.dx = dx
        if(metodo == "adelante"):
            self.met = "adelante"
        elif(metodo == "central"):
            self.met = "central"
        elif(metodo == "extrapolada"):
            self.met = "extrapolada"
        elif(metodo == "segunda"):
            self.met = "segunda"

    def calc(self,x):
        if(self.met == "adelante"):
            return (self.f(x + self.dx) - self.f(x))/self.dx
        elif(self.met == "central"):
            return (self.f(x + (self.dx / 2)) - self.f(x - (self.dx / 2)))/self.dx
        elif(self.met == "extrapolada"):
            f1 = Derivada(self.f, "central", self.dx)
            f2 = Derivada(self.f, "central", (self.dx)/2)
            #print(type(f1), type(f2))
            return (4*f2.calc(x) - f1.calc(x))/3 
        elif(self.met == "segunda"):
            return (self.f(x + self.dx) + self.f(x - self.dx) - 2*self.f(x))/(self.dx*self.dx)

class Zeros:
    def __init__(self, f, metodo, error=1e-4, max_iter=100):
        pass

    def zero(self,vi):
        pass
"""
Adicionalmente el modulo debe tener una sección de ejecución como programa,
donde se muestre un ejemplo del uso de las 2 clases. Esto quiere decir, debe
tener algo del estilo a:
"""
if __name__ == "__main__":
    import math as math
    import numpy as np
    import pylab as plb
    f1 = math.sin
    f2 = math.cos
    f3 = math.exp
    x0 = math.pi
    dx = np.linspace(1,10e-10,1000)
    x = [[],[],[],[]]
    y = [[],[],[],[]]
    j = 0
    while(j<4):
        metodo = "adelante" if j==0 else "central" if j==1 else "extrapolada" if j==2 else "segunda"
        i = 0
        while(i<len(dx)):
            df = Derivada(f1,metodo,dx[i])
            x[j].append(dx[i])
            if(j==3):
                y[j].append(abs(0 - df.calc(x0)))
            else:
                y[j].append(abs(-1 - df.calc(x0)))
            #print("xj = ",x[j])
            #print("yj = ",y[j])
            i += 1
        j += 1

    #line1, = plb.plot(x[0],y[0],linestyle='dashed')
    #line2, = plb.plot(x[1],y[1],linestyle='dashdot')
    #line3, = plb.plot(x[2],y[2],linestyle='dotted')
    line4, = plb.plot(x[3][995:1000],y[3][995:1000])
    #plb.legend((line1, line2, line3, line4),('adelante','central','extrapolada','segunda derivada'))
    #plb.legend((line4,),('segunda',))
    plb.show()
