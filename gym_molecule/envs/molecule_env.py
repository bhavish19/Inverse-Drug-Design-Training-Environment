import gym
from gym import error, spaces, utils
from gym.utils import seeding
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import time
import random
import csv
import itertools
from multiprocessing import Pool
import multiprocessing as mp
#from numba import jit, cuda


class MoleculeEnvironment(gym.Env):
    def __init__(self):
        super().__init__()
        self.currentState = Chem.MolFromSmiles('CCl')
        self.targetState = Chem.MolFromSmiles('CCO')

    def molToString(self, molecule):
        temp = ""
        for a in molecule.GetAtoms():
            temp += a.GetSymbol()
        return temp

    def molToList(self, molecule):
        temp = []
        for a in molecule.GetAtoms():
            temp.append(a.GetSymbol())
        return temp

    def calculateAtomnumbers(self, molecule):
        toReturn = []
        moleculeTemp = self.molToList(molecule)
        symbolList = []
        numberOfAtomsList = []
        k = 0
        while(moleculeTemp):
            moleculeTemp2 = []
            moleculeTemp2 = moleculeTemp.copy()
            moleculeCheck = moleculeTemp[0]
            count = 0
            for i in moleculeTemp:
                if(moleculeCheck == i):
                    count += 1
                    moleculeTemp2.remove(i)
            toReturn.append([])
            toReturn[k].append(moleculeCheck)
            toReturn[k].append(count)
            k += 1
            moleculeTemp = moleculeTemp2
        return toReturn

    def differenceBetweenMols(self, currentMol, targetMol):
        current = self.calculateAtomnumbers(currentMol)
        target = self.calculateAtomnumbers(targetMol)
        symbolList = []
        numberOfAtomsList = []
        eq = []
        for i in target:
            flag = False
            for j in current:
                if (i[0] == j[0]):
                    eq.append(current.index(j))
                    if(i[1] != j[1]):
                        symbolList.append(i[0])
                        temp = i[1] - j[1]
                        numberOfAtomsList.append(temp)
                        flag = True
                        break
            if(flag == False):
                symbolList.append(i[0])
                numberOfAtomsList.append(i[1])
        for i in range(len(current)):
            if i not in eq:
                symbolList.append(current[i][0])
                numberOfAtomsList.append(-current[i][1])
        toReturn = []
        for k in range(len(symbolList)):
            toReturn.append([])
            toReturn[k].append(symbolList[k])
            toReturn[k].append(numberOfAtomsList[k])
        return toReturn

    def combinationsForParallel(self, start, list):
        toReturn = []
        nextMolecule = self.currentStringState + start
        nextValidMolecule = Chem.MolFromSmiles(nextMolecule)
        if nextValidMolecule and not self.currentStringState == nextMolecule:
            toReturn.append(nextMolecule)
        if not list:
            return toReturn
        temp = ""
        count = 0
        for atom, numbers in list:
            if numbers == 1:
                tempList = list.copy()
                del tempList[count]
                start2 = start + atom
                temp2 = self.combinationsForParallel(start2, tempList)
                if temp2:
                    for t in temp2:
                        if t not in toReturn:
                            toReturn.append(t)
            elif numbers > 1:
                tempList = list.copy()
                tempList[count][1] -= 1
                self.combinationsForParallel(start, tempList)
                start2 = start + atom
                temp2 = self.combinationsForParallel(start2, tempList)
                if temp2:
                    for t in temp2:
                        if t not in toReturn:
                            toReturn.append(temp2)
            count += 1
        return toReturn

    def combinations(self, list):
        pool = mp.Pool(mp.cpu_count())
        starters = []
        temp = ""
        starters.append(temp)
        count = 0
        index = 0
        for atom, number in list:
            if 'C' == atom:
                index = count
                for i in range(number):
                    temp += 'C'
                    starters.append(temp)
            count += 1
        tempList = list.copy()
        del tempList[index]
        results = pool.starmap(self.combinationsForParallel, [
                               (row, tempList) for row in starters])
        pool.close()
        return results

    def allValidMolecules(self, atomsList):
        toReturn = []
        atomsPlus = []
        atomsMinus = []
        for i in atomsList:
            if(i[1] > 0):
                atomsPlus.append(i)
        for i in atomsList:
            if(i[1] < 0):
                atomsMinus.append(i)
        p = Pool(4)
        times = range(0, len(atomsPlus)+1)
        values = p.map(self.nLengthCombination(times, atomsPlus))

        p.close()
        return toReturn

    def step(self, nextState, zeroStep):
        self.observations = []
       # file = open('delaney.csv')
       # csv_reader = csv.reader(file, delimiter=",")
       # next(csv_reader, None)
       # rows = list(csv_reader)
       # for i in range(30):
       #     r = random.randint(1, 120)
       #     row = rows[r]
       #     smiles = row[-1]
       #     self.render(smiles)
       #     time.sleep(0)

        #m1 = Chem.MolFromSmiles("C")
        # for i in range (0,10):
        #   smiles = self.randomSmiles(m1)
        #  self.render(smiles)

        modifications = 0  # List for later

        reward = 10  # for now Out Of Scope
        # self.stateReward += reward  # for now Out Of Scope
        #self.currentState = nextState
        self.currentState = Chem.MolFromSmiles('C')
        self.targetState = Chem.MolFromSmiles('OCC(O)CC(O)[C@@H]1C[C@@H]2O[C@@]3(C[C@H](C)[C@@H]2O1)C[C@H](C)[C@@H]4O[C@]%10(C[C@@H]4O3)C[C@H]%11O[C@H]%12[C@H](C)[C@H]%13OC(=O)C[C@H]8CC[C@@H]9O[C@H]7[C@H]6O[C@]5(O[C@H]([C@@H]7O[C@@H]6C5)[C@H]9O8)CC[C@H]%15C/C(=C)[C@H](CC[C@H]%14C[C@@H](C)\C(=C)[C@@H](C[C@@H]%13O[C@H]%12C[C@H]%11O%10)O%14)O%15')
        self.currentStringState = self.molToString(self.currentState)
        atomsList = self.differenceBetweenMols(
            self.currentState, self.targetState)
        validMolecules = self.combinations(atomsList)
        validMolecules = list(itertools.chain.from_iterable(validMolecules))
        validMolecules2 = []
        for i in validMolecules:
            if not isinstance(i, str):
                for j in i:
                    if j not in validMolecules2:
                        validMolecules2.append(j)
            else:
                if i not in validMolecules2:
                    validMolecules2.append(i)
        validMolecules2.sort()
        for i in validMolecules2:
            print(i)
        molecules = 0
        modifications = molecules + reward  # List for later
        return modifications

    def reset(self):
        self.currentState = Chem.MolFromSmiles('C')
        print("Called reset function")

    def render(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        Draw.MolToFile(mol, "molecule.png")
        img = Image.open('molecule.png')
        img.show()
        img.close()
        #img =open('molecule.png','rb').read()

    def seed(self):
        raise NotImplementedError

    def agent(self):
        self.targetState = Chem.MolFromSmiles('nccccc')
        flag = True
        self.stateReward = 0
        nextState = Chem.MolFromSmiles('C')
        modifications = step(self, nextState, True)
        while(flag):
            # Find Next State Here from modifications
            modifications = step(self, nextState, False)

    def getCurrentState(self):
        return self.currentState

    def setCurrentState(self, state):
        self.currentState = state

    def getStateReward(self):
        return self.stateReward

    def setStateReward(self, reward):
        self.stateReward = reward
