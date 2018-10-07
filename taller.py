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
            return (4*f2.calc(x) - f1.calc(x))/3 
        elif(self.met == "segunda"):
            return (self.f(x + self.dx) + self.f(x - self.dx) - 2*self.f(x))/(self.dx*self.dx)

class Zeros:
    def __init__(self, f, metodo, error=1e-4, max_iter=100):
        self.f = f
        self.error = error
        self.max_iter = max_iter
        if(metodo == "newton"):
            self.met = "newton"
        elif(metodo == "bisectriz"):
            self.met = "bisectriz"
        elif(metodo == "interpolacion"):
            self.met = "interpolacion"
        elif(metodo == "newton-sp"):
            self.met = "newton-sp"
        elif(metodo == "fsolve-sp"):
            self.met = "fsolve-sp"
        elif(metodo == "brentq-sp"):
            self.met = "brentq-sp"
        

    def zero(self,vi):
        iteracion = 0
        if(type(vi) == type(3.14159)):
            error = abs(self.f(vi))
            x0 = vi
            df = Derivada(self.f,"extrapolada",0.00000001)

            if(self.met == "newton"):
                while((error >= self.error) and (iteracion <= self.max_iter)):
                    x0 = x0 - (self.f(x0)/df.calc(x0))
                    error = abs(self.f(x0))
                    iteracion += 1

                #return (x0, self.f(x0))
                return x0
                
            elif(self.met == "newton-sp"):
                root = optimize.newton(self.f,vi)
                return root

        elif(type(vi) == type((3,14159))):
            if((self.f(vi[0]) < 0) and (self.f(vi[1]) > 0)):
                    x1 = vi[0]
                    x2 = vi[1]
            elif((self.f(vi[1]) < 0) and (self.f(vi[0]) > 0)):
                    x1 = vi[1]
                    x2 = vi[0]

            if(self.met == "bisectriz"):
                x3 = (x1+x2)/2
                error = abs(self.f(x3))
                iteracion += 1

                while((error >= self.error) and (iteracion <= self.max_iter)):
                    if(self.f(x3) < 0):
                        x1 = x3
                    elif(self.f(x3) > 0):
                        x2 = x3
                    elif(self.f(x3) == 0):
                        return (x3,self.f(x3))
                    x3 = (x1+x2)/2
                    error = abs(self.f(x3))
                    iteracion += 1

                #return (x3,self.f(x3))
                return x3
                
            elif(self.met == "interpolacion"):
                x3 = ((x2*self.f(x1)) - (x1*self.f(x2)))/(self.f(x1)-self.f(x2))
                error = abs(self.f(x3))
                iteracion += 1
                while((error >= self.error) and (iteracion <= self.max_iter)):
                    if(self.f(x3) < 0):
                        x1 = x3
                    elif(self.f(x3) > 0):
                        x2 = x3
                    elif(self.f(x3) == 0):
                        return (x3,self.f(x3))
                    x3 = ((x2*self.f(x1)) - (x1*self.f(x2)))/(self.f(x1)-self.f(x2))
                    error = abs(self.f(x3))
                    iteracion += 1
                #return(x3,self.f(x3))
                return x3
                
            elif(self.met == "fsolve-sp"):
                root = optimize.fsolve(self.f,vi[0])
                return root

            elif(self.met == "brentq-sp"):
                root = optimize.brentq(self.f,vi[0],vi[1])
                return root
        else: 
            print("Método no implementado. \n Los métodos dispinibles son 'newton', 'bisectriz', 'interpolacion', 'newton-sp', 'fsolve-sp', 'brentq-sp'. Donde los terminados en '-sp' significan que son el método usado por el modulo Scipy.")
            return NotImplemented
"""
Adicionalmente el modulo debe tener una sección de ejecución como programa,
donde se muestre un ejemplo del uso de las 2 clases. Esto quiere decir, debe
tener algo del estilo a:
"""
if __name__ == "__main__":
    import math as math
    import numpy as np
    import pylab as plb
    from scipy import optimize
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
            i += 1
        j += 1

    line1, = plb.plot(x[0],y[0],linestyle='dashed')
    line2, = plb.plot(x[1],y[1],linestyle='dashdot')
    line3, = plb.plot(x[2],y[2],linestyle='dotted')
    line4, = plb.plot(x[3],y[3])
    plb.legend((line1, line2, line3, line4),('adelante','central','extrapolada','segunda derivada'))
    plb.show()
    
    x2 = []
    y2 = []
    j = 0
    error = 1e-10
    max_iter = 100
    xi = ((3*math.pi)/4, (5*math.pi)/4)
    while(j<6):
        metodo = "newton" if j==0 else "bisectriz" if j==1 else "interpolacion" if j==2 else "newton-sp" if j==3 else "fsolve-sp" if j==4 else "brentq-sp" 
        raiz = Zeros(f1,metodo,error,max_iter)
        if(j==0 or j==3):
           x1 = raiz.zero(xi[0])
           y1 = j
        else:
            x1 = raiz.zero(xi)
            y1 = j
        x2.append(x1)
        y2.append(y1)
        j += 1

    print(x2)
    print(y2)



