import gym
from gym import error, spaces, utils
from gym.utils import seeding
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import time
import random
import csv


class MoleculeEnvironment(gym.Env):
    def __init__(self):
        super().__init__()

    def step(self):

        file = open('delaney.csv')
        csv_reader = csv.reader(file, delimiter=",")
        next(csv_reader, None)
        rows = list(csv_reader)
        for i in range(30):
            r = random.randint(1, 120)
            row = rows[r]
            smiles = row[-1]
            self.render(smiles)
            time.sleep(0)

       # m1 = Chem.MolFromSmiles("C")
        # for i in range (0,10):
         #   smiles = self.randomSmiles(m1)
          #  self.render(smiles)

    def reset(self):

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
