import gym
import random
import gym_molecule
import os
import time
import sys
import multiprocessing
from rdkit import Chem
from kivy.app import App 
import kivy 
from tkinter import Tk 
from kivy.uix.floatlayout import FloatLayout 
from tkinter.filedialog import askopenfilename
from kivy.uix.boxlayout import BoxLayout
from kivy.config import Config 
from kivy.uix.screenmanager import ScreenManager,Screen
from kivy.lang import Builder
Config.set('graphics', 'resizable', False)  
Config.set('graphics', 'width', '1600')
Config.set('graphics', 'height', '800') 


  
class backgroundTab(Screen):
    
    def spinnerAction(self,text):
        if text == "Upload Action Space":
            self.buttonActionSpace.opacity=1
            self.buttonActionSpace.disabled= False
        else:
            self.buttonActionSpace.opacity=0
            self.buttonActionSpace.disabled= True
    
    def spinnerStartingMolecule(self,text):
        if text == "Random Molecule from List":
            self.buttonStartingMolecule.opacity=1
            self.buttonStartingMolecule.disabled= False
            self.input2.opacity=0
            self.input2.disabled= True
        elif text == "Certain Molecule":
            self.input2.opacity=1
            self.input2.disabled= False
            self.buttonStartingMolecule.opacity=0
            self.buttonStartingMolecule.disabled= True
        else:
            self.buttonStartingMolecule.opacity=0
            self.buttonStartingMolecule.disabled= True
            self.input2.opacity=0
            self.input2.disabled= True

    def actionSpaceFile(self):
        Tk().withdraw()
        filename = askopenfilename()
        print(filename)

    def startingMoleculeFile(self):
        Tk().withdraw()
        filename = askopenfilename()
        print(filename)
    
    def startEnvironment(self):
        pass

    def exit(self):
        exit()

class background(App):  
    def build(self):
        return backgroundTab()
    
    
  

if __name__ == '__main__': 
    background().run()
    multiprocessing.freeze_support()


   
"""
    if __name__ == '__main__':
        root = tk.Tk(className=' TEDD')
    GuiWindow(root).pack()
    root.mainloop()
    multiprocessing.freeze_support()
"""