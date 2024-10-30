import math
from random import gauss
from xml.etree.ElementInclude import include

import zad2
from zad2 import GaussTable


class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def print(self):
        print([self.x, self.y])


class Element:
    def __init__(self, id, jakobian):
        self.id = id
        self.jakobian = jakobian

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
        points = GaussTable(math.sqrt(npc)).returnPoints()
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
        self.printTabs()

    def printTabs(self):
        print("Ksi table: ")
        for i in range(len(self.dN_dksi)):
            print(self.dN_dksi[i])
        print("\nEta table: ")
        for i in range(len(self.dN_deta)):
            print(self.dN_deta[i])


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

def calcH(jakobian, elementUniv, k, npc):
    dN_dx = [[0 for _ in range(4)] for _ in range(npc)]
    dN_dy = [[0 for _ in range(4)] for _ in range(npc)]
    for l, j1 in enumerate(jakobian.J1):
        for x in range(4):
            dN_dx[l][x] = j1[0] * elementUniv.dN_dksi[l][x] +  j1[1] * elementUniv.dN_deta[l][x]
            dN_dy[l][x] = j1[2] * elementUniv.dN_dksi[l][x] +  j1[3] * elementUniv.dN_deta[l][x]
    print("\ndN/dx: ")
    for i in range(len(dN_dx)):
        print(dN_dx[i])
    print("\ndN/dy: ")
    for i in range(len(dN_dy)):
        print(dN_dy[i])

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

    for i in range(len(Hpc)):
        print("Hpc", i+1)
        for j in range(len(Hpc[i])):
            print(Hpc[i][j])
        print()


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

    for i in range(len(H)):
        print(H[i])




def ReadNodesAndElementsFromFile(lines, grid, elementUniv, npc):
    currentLine = 11
    for i in range(currentLine, currentLine + grid.nN):
        nodeData = lines[i].split(", ")
        grid.node.append(Node(float(nodeData[1]), float(nodeData[2])))

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
        grid.element.append(Element(id, jakobian))


file = open('Test1_4_4.txt', 'r')
lines = file.readlines()
file.close()
npcTest = 4
globalData = GlobalData(lines)
elementUniv = ElementUniv(npcTest)
grid = Grid(globalData.nN, globalData.nE)

# ReadNodesAndElementsFromFile(lines, grid, elementUniv, globalData.npc)
# grid.printElementsAndNodes()

# Przyklad z prezentacji
print("Przyklad z prezentacji: ")
gridPrz = Grid(4, 1)
gridPrz.node.append(Node(0.01, -0.01))
gridPrz.node.append(Node(0.025, 0.0))
gridPrz.node.append(Node(0.025, 0.025))
gridPrz.node.append(Node(0.0, 0.025))
jakobian = Jakobian([gridPrz.node[0],
                     gridPrz.node[1],
                     gridPrz.node[2],
                     gridPrz.node[3]],
                    elementUniv, npcTest)

gridPrz.element.append(Element([1, 2, 3, 4], jakobian))
gridPrz.printElementsAndNodes()

calcH(jakobian, elementUniv, 30, npcTest)
# Dla 9 pkt całkowania - z tabelki
