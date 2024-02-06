import csv
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import numpy
from scipy.linalg import lu
from numpy import genfromtxt
from numpy.linalg import linalg

matplotlib.use('TkAgg')

def solver_SOR(A,b):
    tolerancia = 0.00001
    max_iteraciones = 2000
    posicion_valores_090 = [-889, -9, 0, 2, 11, 891]
    posicion_valores_180 = [-3579, -19, 0, 2, 21, 3581]  # No incluyo 1 porque es la diagonal i=j
    posicion_valores_360 = [-17949, -49, 0, 2, 51, 17951]
    return resolver_SOR_optimizado(A, b, w, tolerancia, max_iteraciones, posicion_valores)

def calcularR(X,XAnterior):
    diferencia = [0] * len(X)
    for i in range(len(X)):
        diferencia[i] = X[i] - XAnterior[i]
    maxDif = max([abs(valor) for valor in diferencia])
    return maxDif

def calcularRRelativo(X, XAnterior):
    maxActual = max([abs(valor) for valor in X])
    R = calcularR(X,XAnterior)
    return R / maxActual


def calcularPosicionValoresFila(posicion_valores, i):
    posicion = [x + (i - 1) for x in posicion_valores]
    posicion[:] = [x for x in posicion if (x >= 0 and x < tam_matriz)]
    return posicion

def resolver_GS_original(A, b, tol, max_iteraciones):
    resolver_SOR_original(A, b, 1, tol, max_iteraciones)

def resolver_SOR_original(A, b, w, tolerancia, max_iteraciones):
    start = time.time()
    tam_matriz = len(b)
    X = [3] * tam_matriz  # semilla Tamb

    for iteracion in range(max_iteraciones):
        XAnterior = X.copy()
        for i in range(tam_matriz):
            sum = 0
            for j in range(tam_matriz):
                if j == i:
                    continue
                sum = sum + A[i][j] * X[j]
            X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])

        R = calcularRRelativo(X, XAnterior)
        print("R = " + str(R))
        if R <= tolerancia:
            print("Se llegó a la tolerancia: " + str(tolerancia))
            print("Cantidad de iteraciones: ", iteracion + 1)
            print("R = " + str(R))
            break
    end = time.time()
    print("Tiempo calculo: ", end - start)
    return X


def resolver_SOR_optimizado(A, b, w, tolerancia, max_iteraciones, posicion_valores):
    start = time.time()
    tam_matriz = len(b)
    X = [20] * tam_matriz  # semilla Tamb

    X[0] = b[0]  # La primer fila viene resuelta
    X[tam_matriz - 1] = b[tam_matriz - 1]  # La ultima fila viene resuelta

    for iteracion in range(max_iteraciones):
        XAnterior = X.copy()
        #    for i in range(tam_matriz):
        #        sum = 0
        #        for j in range(tam_matriz):
        #            if j == i:
        #                continue
        #            sum = sum + A[i][j] * X[j]
        #        X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])
        for i in range(1, tam_matriz - 1):
            sum = 0
            posicion_valores_fila = calcularPosicionValoresFila(posicion_valores, i)
            for j in posicion_valores_fila:
                sum = sum + A[i][j] * X[j]
            X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])

        R = calcularRRelativo(X, XAnterior)
        print("R = " + str(R))
        if R <= tolerancia:
            print("Se llegó a la tolerancia: " + str(tolerancia))
            print("Cantidad de iteraciones: ", iteracion + 1)
            print("R = " + str(R))
            break
    end = time.time()
    print("Tiempo calculo: ", end - start)
    return X


def resolver_GS(A, b, tol, max_iteraciones):
    resolver_SOR_optimizado(A, b, 1, tol, max_iteraciones)


def leerCSV(nombreArchivo):
    matriz = list()
    with open(nombreArchivo) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        for fila in csv_reader:
            filaFloat = [float(i) for i in fila]
            filaNumeros = [int(i) for i in filaFloat]
            if len(fila) == 1:
                matriz.append(int(float(fila[0])))
            else:
                matriz.append(filaNumeros)
    return matriz


def cantidadIteracionesSOR(A, b, w, tolerancia, max_iteraciones):
    tam_matriz = len(b)
    X = [20] * tam_matriz  # semilla

    cantidadIteraciones = 0
    for iteracion in range(max_iteraciones):
        XAnterior = X.copy()
        for i in range(tam_matriz):
            sum = 0
            for j in range(tam_matriz):
                if j == i:
                    continue
                sum = sum + A[i][j] * X[j]
            X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])

        R = calcularRRelativo(X, XAnterior)
        if R <= tolerancia:
            cantidadIteraciones = iteracion + 1
            break
    return cantidadIteraciones


def graficarResultadoIteraciones(listaW, resultadosError, iteraciones):
    plt.title("Cantidad iteraciones: " + str(iteraciones))
    plt.xlabel("w")
    plt.ylabel("Cota error")
    plt.plot(listaW, resultadosError, color="green", marker='.')
    plt.show()


def hallarWOptimoIteraciones(A, b, wDesde, wHasta, iteraciones, posicion_valores):
    start = time.time()
    resultadosError = list()
    paso = (wHasta - wDesde) / 10
    listaW = np.arange(wDesde, wHasta, paso)  # start (included), stop (excluded), step
    for w in listaW:
        tam_matriz = len(b)
        X = [20] * tam_matriz  # semilla
        XAnterior = X.copy()

        X[0] = b[0]  # La primer fila viene resuelta
        X[tam_matriz - 1] = b[tam_matriz - 1]  # La ultima fila viene resuelta

        for iteracion in range(iteraciones):
            XAnterior = X.copy()
            for i in range(1, tam_matriz - 1):
                sum = 0
                posicion_valores_fila = calcularPosicionValoresFila(posicion_valores, i)
                for j in posicion_valores_fila:
                    sum = sum + A[i][j] * X[j]
                X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])
        R = calcularRRelativo(X, XAnterior)
        print("w: ", w, " R: ", R)
        resultadosError.append(R)

    end = time.time()
    print("Tiempo calculo: ", end - start)

    graficarResultadoIteraciones(listaW, resultadosError, iteraciones)


def obtenerT(X, numero_filas, numero_columnas):
    T = [[None] * numero_columnas for i in range(numero_filas)]  # Defino la matriz del tamaño requerido

    for i in range(1, ni + 1):
        for j in range(1, nj + 1):
            kx = j + nj * (i - 1)  # fila-columna
            T[i - 1][j - 1] = X[kx - 1]
    return T


def graficarCentroTubo(T, ni, nj):
    plt.title("Temperatura tubo")
    plt.xlabel("ni")
    plt.ylabel("T(ni,nj/2)")

    temperatura = list()
    for i in range(ni):
        temperatura.append(T[i][int(nj / 2)])
    plt.plot(range(ni), temperatura, color="green")
    # plt.savefig("tempTubo.eps", dpi=1200)
    plt.show()


def graficarTuboPolar(T, rext, wt, ni, nj, dr):
    theta = np.linspace(0, 2 * np.pi, ni)
    r = np.linspace(rext - wt, rext, nj)

    R, THETA = np.meshgrid(r, theta)
    Z = np.sin(THETA) * R
    plt.subplot(111, projection='polar')
    plt.pcolormesh(THETA, R, T, cmap='plasma')
    plt.gca().set_rmin(0.0)

    plt.show()


def graficarMatrizT(T, ni, nj):
    plt.rcParams["figure.figsize"] = [1, 4]
    plt.rcParams["figure.autolayout"] = True

    fig, ax = plt.subplots()
    ax.set_title("Matriz T ni=" + str(ni) + "  nj=" + str(nj))
    ax.set_xlabel("ni")
    ax.set_ylabel("nj")
    matrix = T
    ax.matshow(matrix, cmap='twilight')
    # plt.savefig("matrizT.eps", dpi=1200)
    plt.show()


def calcularRadioEspectral(A):
    matriz = np.matrix(A)
    autovalores = linalg.eigvals(A)
    a = list()

    return numpy.abs(autovalores).max()


def hallarTSOR(A1, w):
    A = np.matrix(A1,dtype=np.int32)
    # Descomposicion A = d-l-u
    u = -np.triu(A, 1)  # Separa la parte diagonal superior e invierte los signos
    l = -np.tril(A, -1)  # Separa la parte diagonal inferior e invierte los signos
    d = np.tril(np.triu(A))  # Separa la diagonal
    invertida = numpy.linalg.inv(d - w * l)

    tSOR = np.matmul(invertida,((1 - w) * d + w * u))
    return tSOR

def obtenerErroresSOR(A,B,w,cantidadIteraciones,posicion_valores):
    start = time.time()
    tam_matriz = len(b)
    X = [20] * tam_matriz  # semilla Tamb
    listaErrores = list()
    X[0] = b[0]  # La primer fila viene resuelta
    X[tam_matriz - 1] = b[tam_matriz - 1]  # La ultima fila viene resuelta

    for iteracion in range(cantidadIteraciones):
        XAnterior = X.copy()
        for i in range(1, tam_matriz - 1):
            sum = 0
            posicion_valores_fila = calcularPosicionValoresFila(posicion_valores, i)
            for j in posicion_valores_fila:
                sum = sum + A[i][j] * X[j]
            X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])

        R = calcularRRelativo(X, XAnterior)
        listaErrores.append(R)
        print("R = " + str(R))
        #if R <= tolerancia:
        #    print("Se llegó a la tolerancia: " + str(tolerancia))
        #    print("Cantidad de iteraciones: ", iteracion + 1)
        #    print("R = " + str(R))
        #    break
    end = time.time()
    print("Tiempo calculo: ", end - start)
    return listaErrores

def obtenerErroresGS(A,b,cantidadIteraciones,posicion_valores):
    return obtenerErroresSOR(A,b,1,cantidadIteraciones,posicion_valores)

def graficarErrorRelativoIteraciones(A, b, w, cantidadIteraciones, posicion_valores):
    resultadosErrorSOR = obtenerErroresSOR(A,b,w,cantidadIteraciones,posicion_valores)
    resultadosErrorGS = obtenerErroresGS(A,b,cantidadIteraciones,posicion_valores)
    listaIteraciones = range(1,cantidadIteraciones+1)
    print(listaIteraciones)
    print(len(listaIteraciones))
    print(len(resultadosErrorSOR))
    plt.xlabel("Iteraciones")
    plt.ylabel("Error")
    plt.plot(listaIteraciones, resultadosErrorSOR, color="green", marker='.',label= "SOR")
    plt.plot(listaIteraciones, resultadosErrorGS, color="orange", marker='.', label="GS")
    plt.yscale('log')
    plt.legend()
    plt.show()


def obtenerListasErroresSOR(A,b,w,tolerancia,max_iteraciones,posicion_valores):
    tam_matriz = len(b)
    X = [20] * tam_matriz  # semilla Tamb

    X[0] = b[0]  # La primer fila viene resuelta
    X[tam_matriz - 1] = b[tam_matriz - 1]  # La ultima fila viene resuelta
    listaErrores = list()
    for iteracion in range(max_iteraciones):
        XAnterior = X.copy()
        for i in range(1, tam_matriz - 1):
            sum = 0
            posicion_valores_fila = calcularPosicionValoresFila(posicion_valores, i)
            for j in posicion_valores_fila:
                sum = sum + A[i][j] * X[j]
            X[i] = (1 - w) * X[i] + (b[i] - sum) * (w / A[i][i])

        R = calcularR(X,XAnterior)
        listaErrores.append(R)
        print("R = " + str(R))
        if R <= tolerancia:
            print("Se llegó a la tolerancia: " + str(tolerancia))
            print("Cantidad de iteraciones: ", iteracion + 1)
            print("R = " + str(R))
            break
    return listaErrores

def obtenerListasErroresGrafico(listaErrores):
    if len(listaErrores) % 2 != 0:
        listaErrores.pop()
    listaImparEjeY = listaErrores[1::2]  # Elements from list1 starting from 1 iterating by 2
    listaParEjeX = listaErrores[::2]  # Elements from list1 starting from 0 iterating by 2

    listaImparEjeY[:] = [numpy.log(x) for x in listaImparEjeY]
    listaParEjeX[:] = [numpy.log(x) for x in listaParEjeX]
    return listaParEjeX, listaImparEjeY

def obtenerListasErroresGS(A, b, tolerancia, max_iteraciones, posicion_valores):
    w = 1
    return obtenerListasErroresSOR(A,b,w,tolerancia,max_iteraciones,posicion_valores)

def graficarErrorIteracion(A,b,w,tolerancia,posicion_valores):
    listaErrores = obtenerListasErroresSOR(A,b,w,tolerancia,3000,posicion_valores)
    listaEjeXSOR, listaEjeYSOR = obtenerListasErroresGrafico(listaErrores)

    listaErrores = obtenerListasErroresGS(A, b, tolerancia, 3000, posicion_valores)
    listaEjeXGS, listaEjeYGS = obtenerListasErroresGrafico(listaErrores)

    plt.xlabel("ln(| xk - xk-1|)")
    plt.ylabel("ln(| xk+1 - xk|)")
    plt.plot(listaEjeXSOR, listaEjeYSOR, color="green", marker='.',label= "SOR")
    plt.plot(listaEjeXGS, listaEjeYGS, color="orange", marker='.', label="GS")
    plt.legend()
    plt.show()

def calcularOrdenConvergencia(A, b, w, tolerancia, max_iteraciones, posicion_valores):
    listaErrores = obtenerListasErroresSOR(A, b, w, tolerancia, max_iteraciones, posicion_valores)
    ultimoError = listaErrores.pop()
    anteUltimoError = listaErrores.pop()
    antePenultimoError = listaErrores.pop()
    p = numpy.log(ultimoError/anteUltimoError) / numpy.log(anteUltimoError/antePenultimoError)
    return p

def hallarRadioEspectral(A,w):
    tSOR = hallarTSOR(A, w)

    max_ava = calcularRadioEspectral(tSOR)
    print("Radio Espectral: ", max_ava)

if __name__ == '__main__':
    start = time.time()

    # Inicio ET3

    w = 800
    datos_marea = genfromtxt('marea.csv', delimiter=';')

   resolver_GS_original(A, b, tol, max_iteraciones)

    # Fin ET3


    padron = 100558

    tHot = padron / 100 + 300  # °C
    tAmb = 20
    ni = 180  # 360  # nodos coordenada angular
    nj = 20  # 50  # nodos coordenada radial
    n = ni * nj  # nodos totales

    rExt = 0.250  # radio externo del tubo en metros
    wt = 0.015  # espesor de la pared del tubo en metros
    rInt = rExt - wt  # radio interno del tubo en metros
    dr = wt / (nj - 1)  # delta r de la malla en metros

    #A = leerCSV("A_090_010.csv")
    #b = leerCSV("b_090_010.csv")

    A = leerCSV("A_180_020.csv")
    b = leerCSV("b_180_020.csv")

 #   A = leerCSV("A_360_050.csv")
 #   b = leerCSV("b_360_050.csv")

    posicion_valores_090 = [-889, -9, 0, 2, 11, 891]
    posicion_valores_180 = [-3579, -19, 0, 2, 21, 3581]  # No incluyo 1 porque es la diagonal i=j
    posicion_valores_360 = [-17949, -49, 0, 2, 51, 17951]

    posicion_valores = posicion_valores_180

    tam_matriz = len(b)
    X = [0] * tam_matriz  # semilla
    max_iteraciones = 20000

    tolerancia = 0.00001

    # Resuelvo el SEL A*x=b por el método Gauss-Seidel
    X = resolver_SOR_optimizado(A, b, w, tolerancia, max_iteraciones, posicion_valores)


    # Resuelvo el SEL A*x=b por el método SOR

    while (False):
        opcion = int(input('resolver: 1\n'
                           'Graficar error relativo metodos: 2\n'
                           'Graficar error iteracion: 3\n'
                           'woptimo: 4\n'
                           'Calculo Convergencia p: 5\n'
                           'Calcular Radio Espectral: 6\n'
                           'Input:'))
        if opcion == 1:
            w = float(input('Input w:'))
            X = resolver_SOR_optimizado(A, b, w, tolerancia, max_iteraciones, posicion_valores)
        elif opcion == 2:
            w = float(input('Input w:'))
            cantidadIteraciones = int(input('CantidadIteraciones'))
            graficarErrorRelativoIteraciones(A,b,w,cantidadIteraciones,posicion_valores)
        elif opcion == 3:
            w = float(input('Input w:'))
            tolerancia = 0.00001
            graficarErrorIteracion(A,b,w,tolerancia,posicion_valores)
        elif opcion == 4:
            wDesde = float(input('Input wDesde:'))
            wHasta = float(input('Input wHasta:'))
            cantidadIteraciones = int(input('Input Cantidad iteraciones:'))
         #   except ValueError:
         #       print("Not a number")
            hallarWOptimoIteraciones(A, b, wDesde, wHasta, cantidadIteraciones, posicion_valores)
        elif opcion == 5:
            w = float(input('Input w:'))
            p = calcularOrdenConvergencia(A, b, w, tolerancia, max_iteraciones, posicion_valores)
            print(p)
        elif opcion == 6:
            w = float(input('Input w:'))
            rho = hallarRadioEspectral(A,w)

    w = 1.52  # 1.85
   # X = resolver_SOR(A, b, w, tolerancia, max_iteraciones,posicion_valores_090)
    #  X2 = np.linalg.solve(A, b)
    # Recupero la solución del sistema
    #T = obtenerT(X,ni,nj)
    #  T1 = np.matrix(T)
    #  end = time.time()
    # numpy.savetxt("matrizT.csv", T1, delimiter=";")
    #  print("Tiempo calculo T:",end-start)

    #  graficarCentroTubo(T,ni,nj)

    #  graficarTuboPolar(T,rExt,wt,ni,nj,dr)

    #  graficarMatrizT(T,ni,nj)