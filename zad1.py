class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def print(self):
        print([self.x, self.y])


class Element:
    def __init__(self, id):
        self.id = id

    def print(self):
        print(self.id)


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


def ReadNodesAndElementsFromFile(lines, grid):
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
        grid.element.append(Element(id))


file = open('Test1_4_4.txt', 'r')
lines = file.readlines()
file.close()
globalData = GlobalData(lines)
grid = Grid(globalData.nN, globalData.nE)

ReadNodesAndElementsFromFile(lines, grid)
grid.printElementsAndNodes()
