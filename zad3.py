import math
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
        print(self.jakobian.J[0])
        print(self.jakobian.J[1])
        print("Odwrocony Jakobian: ")
        print(self.jakobian.J1[0])
        print(self.jakobian.J1[1])
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
        point = 1 / math.sqrt(3)
        points = [
            # ksi, eta
            -point, -point,
            point, -point,
            -point, point,
            point, point
        ]
        for i in range(0, npc):
            self.dN_dksi.append([-0.25 * (1 - points[2 * i + 1]),
                                 0.25 * (1 - points[2 * i + 1]),
                                 0.25 * (1 + points[2 * i + 1]),
                                 -0.25 * (1 + points[2 * i + 1])])

            self.dN_deta.append([-0.25 * (1 - points[2 * i]),
                                 -0.25 * (1 + points[2 * i]),
                                 0.25 * (1 + points[2 * i]),
                                 0.25 * (1 - points[2 * i])])
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
        self.J = [[0, 0], [0, 0]]
        self.J1 = [[], []]
        for i in range(npc):
            self.J[0][0] += elementUniv.dN_dksi[0][i] * nodeEl[i].x
            self.J[0][1] += elementUniv.dN_dksi[0][i] * nodeEl[i].y
            self.J[1][0] += elementUniv.dN_deta[0][i] * nodeEl[i].x
            self.J[1][1] += elementUniv.dN_deta[0][i] * nodeEl[i].y
        for i in range(2):
            if abs(self.J[i][0]) < 0.00001:
                self.J[i][0] = 0
            if abs(self.J[i][1]) < 0.00001:
                self.J[i][1] = 0
        self.det = self.J[0][0] * self.J[1][1] - (self.J[0][1] * self.J[1][0])

        self.J1[0].append((1 / self.det) * self.J[1][1])
        self.J1[0].append((1 / self.det) * self.J[0][1])
        self.J1[1].append((1 / self.det) * self.J[1][0])
        self.J1[1].append((1 / self.det) * self.J[0][0])

def calcH(jakobian, elementUniv, k, npc):
    dN_dx = [[0 for _ in range(npc)] for _ in range(npc)]
    dN_dy = [[0 for _ in range(npc)] for _ in range(npc)]
    for y in range(npc):
        for x in range(npc):
            dN_dx[y][x] = jakobian.J1[0][0] * elementUniv.dN_dksi[y][x] +  jakobian.J1[0][1] * elementUniv.dN_deta[y][x]
            dN_dy[y][x] = jakobian.J1[1][0] * elementUniv.dN_dksi[y][x] +  jakobian.J1[1][1] * elementUniv.dN_deta[y][x]
    print("\ndN/dx: ")
    for i in range(len(dN_dx)):
        print(dN_dx[i])
    print("\ndN/dy: ")
    for i in range(len(dN_dy)):
        print(dN_dy[i])

#     Mnożenie transponownaych i zwyklych
    HpcX = [[[0 for _ in range(npc)] for _ in range(npc)] for _ in range(npc)]
    HpcY= [[[0 for _ in range(npc)] for _ in range(npc)] for _ in range(npc)]
    for integrationPoint in range(npc):
        for x in range(npc):
            for y in range(npc):
                HpcX[integrationPoint][x][y] = dN_dx[integrationPoint][x] * dN_dx[integrationPoint][y]
                HpcY[integrationPoint][x][y] = dN_dy[integrationPoint][x] * dN_dy[integrationPoint][y]

# Obliczanie Hpc dla każdego punktu całkowania

    Hpc = [[[0 for _ in range(npc)] for _ in range(npc)] for _ in range(npc)]
    for integrationPoint in range(npc):
        for y in range(npc):
            for x in range(npc):
                Hpc[integrationPoint][y][x] =k*(HpcX[integrationPoint][y][x] + HpcY[integrationPoint][y][x]) * 0.000156

    for i in range(len(Hpc)):
        print("Hpc", i+1)
        for j in range(len(Hpc[i])):
            print(Hpc[i][j])
        print()


# Obliczenie H
    gauss = GaussTable(math.sqrt(npc))
    weights = gauss.returnWeights()
    H = [[0 for _ in range(npc)] for _ in range(npc)]
    for y in range(npc):
        for x in range(npc):
            for w in range(npc):
                H[x][y] += Hpc[w][y][x] * weights[int(w // math.sqrt(npc))]

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
globalData = GlobalData(lines)
elementUniv = ElementUniv(globalData.npc)
grid = Grid(globalData.nN, globalData.nE)

# ReadNodesAndElementsFromFile(lines, grid, elementUniv, globalData.npc)
# grid.printElementsAndNodes()

# Przyklad z prezentacji
print("Przyklad z prezentacji: ")
gridPrz = Grid(4, 1)
gridPrz.node.append(Node(0.0, 0.0))
gridPrz.node.append(Node(0.025, 0.0))
gridPrz.node.append(Node(0.025, 0.025))
gridPrz.node.append(Node(0.0, 0.025))
jakobian = Jakobian([gridPrz.node[0],
                     gridPrz.node[1],
                     gridPrz.node[2],
                     gridPrz.node[3]],
                    elementUniv, 4)

gridPrz.element.append(Element([1, 2, 3, 4], jakobian))
gridPrz.printElementsAndNodes()

calcH(jakobian, elementUniv, 30, 4)
# Dla 9 pkt całkowania - z tabelki
