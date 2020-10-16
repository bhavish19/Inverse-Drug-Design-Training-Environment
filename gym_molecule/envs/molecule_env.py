#Importing Libraries
import gym
from gym import error, spaces, utils
from gym.utils import seeding
from rdkit import Chem
import networkx as nx
import random
from multiprocessing import Pool
import multiprocessing as mp
from rdkit.Chem import rdChemReactions
import time
from rdkit.Chem import Draw
from tkinter import Tk 
from tkinter import messagebox

# Creating global reward
reward = float(0)


class MoleculeEnvironment(gym.Env):

    # Initilazer
    def __init__(self):
        super().__init__()


    def init(self, startingMoleculeList, actionSpaceFile):
        """
        Initilazing The Environment

        Args:
            startingMoleculeList (list): List of Starting Molecules
            actionSpaceFile (str): Action Space File Path
        """
        # For parallelThreadization and GUI
        mp.freeze_support()
        self.startingMoleculeList = startingMoleculeList
        # Choosing Random Molecule From List and Making it Current State to define currentState
        self.currentState = Chem.MolFromSmiles(
            random.choice(startingMoleculeList))
        # Reading Action Space File and Creating List of it
        self.reactionsList = []     
  
        with open(actionSpaceFile) as f:
            mylist = [line.rstrip('\n') for line in f]
        try:
            for line in mylist:
                reaction = rdChemReactions.ReactionFromSmarts(line.strip())
                self.reactionsList.append(reaction)
        except:
            print("Error in one of the Reactions")



    def step(self, action):
        """
        Step Method that takes steps (Makes Changes) in environment

        Args:
            action (networkx.Graph()): Action That Will be Applied to currentState

        Returns:
            tuple: observations, reward, Done; Returns observations for next State, Reward for Agents Chosen Action, 
            Done is True if No more Modifications can be made
        """
        validModifications = []
        
        #Changes Current State
        self.take_action(action)
        Done = False

        #Iterates Through reactionList(ActionSpace) and Checks if any Reactants Match to currentState
        #If they Match Checks Products if they are Chemically Valid
        #If they Are Chemically valid adds them to validModifications (list) 
        for reactions in self.reactionsList:
            for reactant in range(reactions.GetNumReactantTemplates()):
                if Chem.MolToSmiles(reactions.GetReactantTemplate(reactant)) == Chem.MolToSmiles(self.currentState):
                    for i in range(reactions.GetNumProductTemplates()):
                        product = reactions.GetProductTemplate(i)
                        if product:
                            validModifications.append(product)

        #Environment can't make any modifications to current State"
        if(len(validModifications) == 0):            
            Done = True

        #Works only when used without GUI
        #Created Molecular Graphs of Molecules in parallelThread
        #p = mp.Pool(mp.cpu_count())
        #observations = p.map(self.moltoGraph, validModifications)
        #p.close

        #Converts Valid Modifications from validModifications (list) to Moleculer Graphs and adds them to observations (list)
        observations = []
        for i in validModifications:
            observations.append(self.moltoGraph(i))

        #Creating tuple to return
        tuple = (observations, self.reward(), Done)
        return tuple


    def reset(self):
        """
        Reseting Environment to the Beginning

        Returns:
            [networkx.Graph()]: currentState in Molecule Graph Format
        """
        global reward
        #Reset Reward
        reward = 0
        #Select Molecule Random from startingMoleculeList (list) and make it currentState
        tempMol= Chem.MolFromSmiles(random.choice(self.startingMoleculeList))
        count = 0
        while not tempMol or count > 20:
            tempMol= Chem.MolFromSmiles(random.choice(self.startingMoleculeList))
            count+=1
        if count>20:
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Please upload Proper Action List")
            tempWindow.destroy()
            exit()

        self.currentState = tempMol
        return self.moltoGraph(self.currentState)


    def render(self):
        """
        Renders CurrentState and saves it as RenderingImage.png to System Folder for GUI to catch
        """
        try:
            Chem.Draw.MolToFile(self.currentState, "RenderingImage.png")
        except:
            pass


    def seed(self):
        """
        FOR USER TO IMPLEMENT
        Not Implemented

        Raises:
            NotImplementedError:
        """
        raise NotImplementedError

 
    def reward(self):
        """
        Dummy Rewarding System
        THIS METHOD IS FOR USER TO IMPLEMENT REWARDING SYSTEM

        Returns:
            [float]: Reward
        """        
        global reward
        #Adds 10 points no matter what action is chosen by Agent
        reward +=10
        return reward


    def take_action(self, action):
        """Takes action to Modify currentState

        Args:
            action (networkx.Graph()): Action That Will be Applied to currentState
        """        
        self.currentState = self.graphToMol(action)

    """
    Out Of Scope
    Converting Molecules to Molecular Graph using networkx library
    For Methods moltoGraph(mol), graphToMol(molGraph)
    Please Check the Referece Github For This Code
    https://github.com/dakoner/keras-molecules/tree/dbbb790e74e406faa70b13e8be8104d9e938eba2
    Might Be Modified to corporate the methods with this Project
    When Creating the Graph or Creating Molecule from Graph
    Doesnt Take into Account Protonation State, ...
    Also Some Bond Types Triple, Dative, ...
    This methods can easily be extended if needed
    |  |  |  |  |  |  |  |   |  |  |  |  |  |  |  |
    \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/
    """
    def moltoGraph(self, mol):
        """
        Converts Molecule to Moleculer Graph

        Args:
            mol (Chem.mol): Molecule in Rdkit chem library

        Returns:
            networkx.Graph(): Molecular Graph that has been converted from Molecule
        """
        molGraph = nx.Graph()

        for atom in mol.GetAtoms():
            molGraph.add_node(atom.GetIdx(),
                              atomic_num=atom.GetAtomicNum(),
                              formal_charge=atom.GetFormalCharge(),
                              chiral_tag=atom.GetChiralTag(),
                              hybridization=atom.GetHybridization(),
                              num_explicit_hs=atom.GetNumExplicitHs(),
                              is_aromatic=atom.GetIsAromatic())
        for bond in mol.GetBonds():
            molGraph.add_edge(bond.GetBeginAtomIdx(),
                              bond.GetEndAtomIdx(),
                              bond_type=bond.GetBondType())
        return molGraph

    def graphToMol(self, molGraph):
        """
        Converts Moleculer Graph to Molecule

        Args:
            molGraph (networkx.Graph()): Molecular Graph

        Returns:
            [Chem.mol]: Molecule in Rdkit chem library that is converted back from molecular graph
        """
        mol = Chem.RWMol()
        atomic_nums = nx.get_node_attributes(molGraph, 'atomic_num')
        chiral_tags = nx.get_node_attributes(molGraph, 'chiral_tag')
        formal_charges = nx.get_node_attributes(molGraph, 'formal_charge')
        node_is_aromatics = nx.get_node_attributes(molGraph, 'is_aromatic')
        node_hybridizations = nx.get_node_attributes(molGraph, 'hybridization')
        num_explicit_hss = nx.get_node_attributes(molGraph, 'num_explicit_hs')
        node_to_idx = {}
        for node in molGraph.nodes():
            a = Chem.Atom(atomic_nums[node])
            a.SetChiralTag(chiral_tags[node])
            a.SetFormalCharge(formal_charges[node])
            a.SetIsAromatic(node_is_aromatics[node])
            a.SetHybridization(node_hybridizations[node])
            a.SetNumExplicitHs(num_explicit_hss[node])
            idx = mol.AddAtom(a)
            node_to_idx[node] = idx

        bond_types = nx.get_edge_attributes(molGraph, 'bond_type')
        for edge in molGraph.edges():
            first, second = edge
            ifirst = node_to_idx[first]
            isecond = node_to_idx[second]
            bond_type = bond_types[first, second]
            mol.AddBond(ifirst, isecond, bond_type)
        try:
            Chem.SanitizeMol(mol)
        except:
            pass

        return mol