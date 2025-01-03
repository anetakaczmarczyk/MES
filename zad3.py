import itertools
import math
from random import gauss
from xml.etree.ElementInclude import include

import zad2
from zad2 import GaussTable


class Node:
    def __init__(self, x, y, isHbc):
        self.x = x
        self.y = y
        self.isHbc = isHbc

    def print(self):
        print([self.x, self.y])


class Element:
    def __init__(self, id, jakobian, node):
        self.id = id
        self.jakobian = jakobian
        self.node = node
        self.matrixH = []
        self.matrixC = []
        self.matrixHBC = []
        self.P=[]

    def print(self):
        print("\nnr id: ", self.id)
        print("Jakobian: ")
        print(self.jakobian.J)
        print("Odwrocony Jakobian: ")
        print(self.jakobian.J1)
        print("Det: ", self.jakobian.det)


class Grid:
    def __init__(self, nN, nE):
        self.nN = nN
        self.nE = nE
        self.node = []
        self.element = []

    def printElementsAndNodes(self):
        print("\nNodes: ")
        for i in range(self.nN):
            self.node[i].print()
        print("\nElements: ")
        for i in range(self.nE):
            self.element[i].print()


class GlobalData:

    def __init__(self, lines):
        self.SimulationTime = int(lines[0].split(" ")[1])
        self.SimulationStepTime = int(lines[1].split(" ")[1])
        self.Conductivity = int(lines[2].split(" ")[1])
        self.Alfa = int(lines[3].split(" ")[1])
        self.Tot = int(lines[4].split(" ")[1])
        self.InitialTemp = int(lines[5].split(" ")[1])
        self.Density = int(lines[6].split(" ")[1])
        self.SpecificHeat = int(lines[7].split(" ")[1])
        self.nN = int(lines[8].split(" ")[2])
        self.nE = int(lines[9].split(" ")[2])
        self.npc = 4


class ElementUniv:
    def __init__(self, npc):
        self.dN_dksi = []
        self.dN_deta = []
        self.Npc = []
        points = GaussTable(math.sqrt(npc)).returnPoints()
        self.surface = [Surface() for _ in range(4)]
        for ksi in points:
            for eta in points:
                self.dN_dksi.append([-0.25 * (1 - eta),
                                     0.25 * (1 - eta),
                                     0.25 * (1 + eta),
                                     -0.25 * (1 + eta)])

                self.dN_deta.append([-0.25 * (1 - ksi),
                                     -0.25 * (1 + ksi),
                                     0.25 * (1 + ksi),
                                     0.25 * (1 - ksi)])
                self.Npc.append([0.25 * (1 - ksi) * (1 - eta),
                 0.25 * (1 + ksi) * (1 - eta),
                 0.25 * (1 + ksi) * (1 + eta),
                 0.25 * (1 - ksi) * (1 + eta)])

        for point in points:
            self.surface[0].addN(point,-1)  #down
            self.surface[1].addN(1,point)  #right
            self.surface[2].addN(point,1)  #up
            self.surface[3].addN(-1,point)  #left




    def printTabs(self):
        print("Ksi table: ")
        for i in range(len(self.dN_dksi)):
            print(self.dN_dksi[i])
        print("\nEta table: ")
        for i in range(len(self.dN_deta)):
            print(self.dN_deta[i])

class Surface:
    def __init__(self):
        self.N = []

    def addN(self, ksi, eta):
        list = [0.25 * (1 - ksi) * (1 - eta),
                0.25 * (1 + ksi) * (1 - eta),
                0.25 * (1 + ksi) * (1 + eta),
                0.25 * (1 - ksi) * (1 + eta)]
        self.N.append(list)


class Jakobian:

    def __init__(self, nodeEl, elementUniv, npc):
        self.J =[[0.0 for _ in range(4)] for _ in range(npc)]
        self.J1 = [[0.0 for _ in range(4)] for _ in range(npc)]
        self.det = []
        for l, j in enumerate(self.J):
            for i in range(4):
                j[0] += elementUniv.dN_dksi[l][i] * nodeEl[i].x
                j[1] += elementUniv.dN_dksi[l][i] * nodeEl[i].y
                j[2] += elementUniv.dN_deta[l][i] * nodeEl[i].x
                j[3] += elementUniv.dN_deta[l][i] * nodeEl[i].y
            for i in range(4):
                if abs(j[i]) < 0.00001:
                    j[i] = 0

            self.det.append(j[0] * j[3] - (j[1] * j[2]))
            self.J1[l][0] = (1 / self.det[l]) * j[3]
            self.J1[l][1] = (1 / self.det[l]) * -j[1]
            self.J1[l][2] = (1 / self.det[l]) * -j[2]
            self.J1[l][3] = (1 / self.det[l]) * j[0]

def gauss_elimination(A, b):
    n = len(A)
    for i in range(n):
        A[i].append(b[i])

    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[max_row] = A[max_row], A[i]

        diag = A[i][i]
        for k in range(i, n + 1):
            A[i][k] /= diag
        for j in range(i + 1, n):
            factor = A[j][i]
            for k in range(i, n + 1):
                A[j][k] -= factor * A[i][k]
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = A[i][n] - sum(A[i][j] * x[j] for j in range(i + 1, n))
    return x


def calcHbc(surface, nodes, npc, alpha, temperature):
    weights = GaussTable(math.sqrt(npc)).returnWeights()
    Hbc = [[0.0 for _ in range(4)] for _ in range(4)]
    P = [0.0 for _ in range(4)]

    for i in range(4):
        if nodes[i].isHbc and nodes[(i + 1) % 4].isHbc:
            length = math.sqrt((nodes[(i + 1) % 4].x - nodes[i].x) ** 2 + (nodes[(i + 1) % 4].y - nodes[i].y) ** 2)
            detJ = length / 2
            for n in range(int(math.sqrt(npc))):
                for j in range(4):
                    for k in range(4):
                        Hbc[j][k] += detJ * alpha * (surface[i].N[n][j] * surface[i].N[n][k]) * weights[n]
                    P[j] += alpha * weights[n] * surface[i].N[n][j] * temperature *detJ
    return Hbc, P



def calcH(el, jakobian, elementUniv, k, npc, surface, nodes, alpha):
    dN_dx = [[0 for _ in range(4)] for _ in range(npc)]
    dN_dy = [[0 for _ in range(4)] for _ in range(npc)]
    for l, j1 in enumerate(jakobian.J1):
        for x in range(4):
            dN_dx[l][x] = j1[0] * elementUniv.dN_dksi[l][x] +  j1[1] * elementUniv.dN_deta[l][x]
            dN_dy[l][x] = j1[2] * elementUniv.dN_dksi[l][x] +  j1[3] * elementUniv.dN_deta[l][x]

#     Mnożenie transponownaych i zwyklych
    HpcX = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    HpcY= [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    for integrationPoint in range(npc):
        for x in range(4):
            for y in range(4):
                HpcX[integrationPoint][x][y] = dN_dx[integrationPoint][x] * dN_dx[integrationPoint][y]
                HpcY[integrationPoint][x][y] = dN_dy[integrationPoint][x] * dN_dy[integrationPoint][y]


# Obliczanie Hpc dla każdego punktu całkowania

    Hpc = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(npc)]
    for integrationPoint in range(npc):
        for y in range(4):
            for x in range(4):
                Hpc[integrationPoint][y][x] =k*(HpcX[integrationPoint][y][x] + HpcY[integrationPoint][y][x]) * jakobian.det[integrationPoint]



# Obliczenie H
    gauss = GaussTable(math.sqrt(npc))
    weights = gauss.returnWeights()
    H = [[0 for _ in range(4)] for _ in range(4)]
    for i in range(len(Hpc)):
        w1 = i % int(math.sqrt(npc))
        w2 = (i//int(math.sqrt(npc))) % int(math.sqrt(npc))
        for x in range(4):
            for y in range(4):
                H[x][y] +=  Hpc[i][x][y] * weights[w1] * weights[w2]
    el.matrixHBC, el.P = calcHbc(surface, nodes, npc, alpha, globalData.Tot)
    for x in range(4):
        for y in range(4):
            H[x][y] += el.matrixHBC[x][y]

    # for i in range(len(H)):
    #     print(H[i])
    # print()
    return H

def calcC(el, elementUniv, jakobian,  npc, c, ro):
    sum = [[0.0 for _ in range(4)] for _ in range(4)]
    weights = GaussTable(math.sqrt(npc)).returnWeights()
    combinations = list(itertools.product(weights, repeat=2))
    for i in range(npc):
        for x in range(4):
            for y in range(4):
                sum[x][y] += elementUniv.Npc[i][x] * elementUniv.Npc[i][y] * c * ro * jakobian.det[i] * combinations[i][0] * combinations[i][1]
    return sum
class GlobalSystemOfEquations :
    def __init__(self, numberOfNodes):
        self.numberOfNodes = numberOfNodes
        self.globalMartixH = [[0 for _ in range(numberOfNodes)] for _ in range(numberOfNodes)]
        self.globalMartixC = [[0 for _ in range(numberOfNodes)] for _ in range(numberOfNodes)]
        self.globalVectorP = [0 for _ in range(numberOfNodes)]

    def calcTemps(self, time, step, initialTemp):

        currTemp = [initialTemp for _ in range(self.numberOfNodes)]


        for currTime in range(step, time+1, step):
            start = [[0.0 for _ in range(self.numberOfNodes)] for _ in range(self.numberOfNodes)]
            end = []
            C = [[0.0 for _ in range(self.numberOfNodes)] for _ in range(self.numberOfNodes)]
            for i in range(self.numberOfNodes):
                for j in range(self.numberOfNodes):
                    start[i][j] = self.globalMartixH[i][j] + self.globalMartixC[i][j] / step
                    C[i][j] = self.globalMartixC[i][j]/step
            for i in range(self.numberOfNodes):
                row = C[i]
                result = 0
                for j in range(self.numberOfNodes):
                    result += row[j] * currTemp[j]
                end.append(result + self.globalVectorP[i])

            print(currTime)
            x = gauss_elimination(start, end)
            print(min(x), max(x))
            currTemp = x



def calcGlobalH(globalMartixH, allElements):
    for element in allElements:
        for i in range(len(element.matrixH)):
            for j in range(4):
                globalMartixH[element.id[i]-1][element.id[j]-1] += element.matrixH[i][j]

def calcGlobalC(globalMartixC, allElements):
    for element in allElements:
        for i in range(len(element.matrixC)):
            for j in range(4):
                globalMartixC[element.id[i] - 1][element.id[j] - 1] += element.matrixC[i][j]

    # for i in range(len(globalMartixC)):
    #     print(globalMartixC[i])

def calcGlobalP(globalVectorP, allElements):
    for element in allElements:
        for i in range(len(element.P)):
            globalVectorP[element.id[i]-1] += element.P[i]

    # for i in range(len(globalVectorP)):
    #     print(globalVectorP[i])

def ReadNodesAndElementsFromFile(lines, grid, elementUniv, npc):
    currentLine = 11
    BCNodes = lines[currentLine + grid.nN +1 + grid.nE + 1].split(", ")
    BCNodes = [item.strip() for item in BCNodes]
    for i in range(currentLine, currentLine + grid.nN):
        nodeData = lines[i].split(", ")
        grid.node.append(Node(float(nodeData[1]), float(nodeData[2]),nodeData[0].strip() in BCNodes))
    currentLine += globalData.nN + 1

    for i in range(currentLine, currentLine + grid.nE):
        elementData = lines[i].split(", ")
        id = []
        for j in range(1, 5):
            id.append(int(elementData[j]))
        jakobian = Jakobian([grid.node[id[0] - 1],
                             grid.node[id[1] - 1],
                             grid.node[id[2] - 1],
                             grid.node[id[3] - 1]],
                            elementUniv, npc)
        node = [grid.node[id[0] - 1],
                grid.node[id[1] - 1],
                grid.node[id[2] - 1],
                grid.node[id[3] - 1]]
        grid.element.append(Element(id, jakobian, node))


# file = open('test.txt', 'r')
# file = open('Test1_4_4.txt', 'r')
file = open('Test2_4_4_MixGrid.txt', 'r')
# file = open('Test3_31_31_kwadrat.txt', 'r')
lines = file.readlines()
file.close()
globalData = GlobalData(lines)
elementUniv = ElementUniv(globalData.npc)
grid = Grid(globalData.nN, globalData.nE)

ReadNodesAndElementsFromFile(lines, grid, elementUniv, globalData.npc)
# grid.printElementsAndNodes()
for el in grid.element:
    el.matrixH = calcH(el, el.jakobian, elementUniv, 25, globalData.npc, elementUniv.surface, el.node, globalData.Alfa)
    el.matrixC = calcC(el, elementUniv, el.jakobian, globalData.npc, globalData.SpecificHeat, globalData.Density)
globalSystemOfEquations = GlobalSystemOfEquations(globalData.nN)
calcGlobalH(globalSystemOfEquations.globalMartixH, grid.element)
calcGlobalP(globalSystemOfEquations.globalVectorP, grid.element)
calcGlobalC(globalSystemOfEquations.globalMartixC, grid.element)

globalSystemOfEquations.calcTemps(globalData.SimulationTime, globalData.SimulationStepTime, globalData.InitialTemp)
# x = gauss_elimination(globalSystemOfEquations.globalMartixH, globalSystemOfEquations.globalVectorP)
# for i in range(len(x)):
#     print(x[i])
# Przyklad z prezentacji
# print("Przyklad z prezentacji: ")
# gridPrz = Grid(4, 1)
# gridPrz.node.append(Node(0.01, -0.01))
# gridPrz.node.append(Node(0.025, 0.0))
# gridPrz.node.append(Node(0.025, 0.025))
# gridPrz.node.append(Node(0.0, 0.025))
# jakobian = Jakobian([gridPrz.node[0],
#                      gridPrz.node[1],
#                      gridPrz.node[2],
#                      gridPrz.node[3]],
#                     elementUniv, npcTest)
#
# gridPrz.element.append(Element([1, 2, 3, 4], jakobian))
# gridPrz.printElementsAndNodes()
#
# calcH(jakobian, elementUniv, 30, npcTest)
# Dla 9 pkt całkowania - z tabelki
