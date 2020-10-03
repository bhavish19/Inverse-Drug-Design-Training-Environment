import gym
from gym import error, spaces, utils
from gym.utils import seeding
from rdkit import Chem
from itertools import chain, combinations


class MoleculeEnvironment(gym.Env):
    def __init__(self):
        super().__init__()
        self.currentState = Chem.MolFromSmiles('C')
        self.targetState = Chem.MolFromSmiles('cchhhh')

    def step(self, nextState, zeroStep):
        modifications = 0  # List for later
        if(zeroStep == False):
            self.currentState = nextState
        reward = 10  # for now Out Of Scope
        self.stateReward += reward  # for now Out Of Scope

        atomsList = differenceBetweenMols(self, self.currentState, self.targetState)
        validMolecules = allValidMolecules(self, atomsList)
        modifications = modificationsToCurrentState(self,validMolecules)
        #molecules = 0
        #modifications = molecules + reward  # List for later
        return modifications

    def reset(self):
        self.currentState = Chem.MolFromSmiles('C')

    def render(self):
        raise NotImplementedError

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

    def molToString(self, molecule):
        temp = ""
        moleculeWithHydrogen = Chem.AddHs(molecule)
        for a in moleculeWithHydrogen.GetAtoms():
            temp += a.GetSymbol()
        return temp

    def modificationsToCurrentState(self, validMolecules):
        current = molToString(self, self.currentState)
        toReturn = []
        for i in validMolecules:
            if i[1]== 1:
                temp = current + i[0]
                mol = Chem.MolFromSmiles(temp)
                toReturn.append(mol)
            

        return toReturn

    def calculateAtomnumbers(self, molecule):
        moleculeTemp = molToString(self, molecule)
        symbolList = []
        numberOfAtomsList = []
        for i in moleculeTemp:
            count = 0
            for j in symbolList:
                if(i == j):
                    numberOfAtomsList[count] += 1
                    break
                symbolList.append(i)
                numberOfAtomsList.append(1)
                count += 1
        toReturn = []

        for c in len(symbolList):
            toReturn[c] = []
            toReturn[c][0] = symbolList[c]
            toReturn[c][1] = numberOfAtomsList[c]

        return toReturn

    def differenceBetweenMols(self, currentMol, targetMol):
        current = calculateAtomnumbers(self, currentMol)
        target = calculateAtomnumbers(self, targetMol)
        symbolList = []
        numberOfAtomsList = []
        for i in target:
            for j in current:
                if (i[0] == j[0]):
                    if(i[1] != j[1]):
                        symbolList.append(i[0])
                        temp = i[1] - j[1]
                        numberOfAtomsList.append(temp)
                        break
        toReturn = []
        for c in len(symbolList):
            toReturn[c] = []
            toReturn[c][0] = symbolList[c]
            toReturn[c][1] = numberOfAtomsList[c]
        return toReturn

    def allMolecules(self, atomsList):
        toReturn = []
        atomsPlus = []
        atomsMinus = []
        for i in atomsList:
            if(i[1] > 0):
                for j in i[1]:
                    atomsPlus.append(i[0])
        for i in atomsList:
            if(i[1] < 0):
                for j in i[1]:
                    atomsMinus.append(i[0])

        moleculesPlus = allCombinations(self, atomsPlus)
        moleculesMinus = allCombinations(self, atomsMinus)
        toReturn[0] = moleculesPlus
        toReturn[1] = moleculesMinus
        return toReturn

    def allValidMolecules(self, atomsList):
        moleculesList = allMolecules(self, atomsList)
        toReturn = []
        c = 0
        for i in moleculesList[0]:
            if(checkChemicalValidity(self, i)):
                toReturn[c] = []
                toReturn[c][0] = i
                toReturn[c][1] = 1
                c += 1
        for i in moleculesList[1]:
            if(checkChemicalValidity(self, i)):
                toReturn[c] = []
                toReturn[c][0] = i
                toReturn[c][1] = -1
                c += 1

        return toReturn

    def allCombinations(self, atomList):
        return chain(*map(lambda x: combinations(atomList, x), range(0, len(atomList)+1)))

    def checkChemicalValidity(self, molecule):
        # implicitly performs sanitization
        molecule2 = Chem.MolFromSmiles(molecule)
        if molecule2:
            return True
        else:
            return False

    def getCurrentState(self):
        return self.currentState

    def setCurrentState(self, state):
        self.currentState = state

    def getStateReward(self):
        return self.stateReward

    def setStateReward(self, reward):
        self.stateReward = reward
