import gym
import random
import gym_molecule
import tkinter as tk
from tkinter import Text, filedialog, ttk, PhotoImage, Label
import os
import time
import sys
from rdkit import Chem


class GuiWindow(tk.Frame):

    def __init__(self, parent):
        self.env = gym.make("gym_molecule:molecule-v0")
        tk.Frame.__init__(self, parent)
        self.photo = PhotoImage(file="demoMolecule.png")
        self.photo = self.photo.subsample(8, 8)
        parent.geometry('500x300')
        self.frame = tk.Frame(root, bg="white")
        self.frame.pack(expand='true', fill="both")
        self.quitButton = tk.Button(self.frame, text="Quit", command=sys.exit)
        self.quitButton.pack(side='left')
        self.restartButton = tk.Button(
            self.frame, text="Restart", command=self.restart_program)
        self.restartButton.pack(side='right')
        self.T = Text(self.frame, height=3)
        Welcome = """Welcome! \nPlease upload a target molecule\nspecified as SMILES"""
        self.T.configure(font=("Courier", 16,))
        self.T.pack()
        self.T.insert(tk.END, Welcome)
        self.addFileButton = tk.Button(
            self.frame, text="Add File", fg="black", bg='blue', command=self.addFile)
        self.addFileButton.pack()

    def addFile(self):
        filename = filedialog.askopenfilename(initialdir="/Users/constantnel/Downloads",
                                              title="Select file", filetypes=(("txt files", "*.txt"), ("all files", "*.*")))
        if(os.path.getsize(filename) > 0):
            file = open(filename, "r")
            smiles = file.readline()
            self.env.render(smiles)
            self.buildNewButtonsAndText()
        else:
            self.buildErrorState()

    def buildNewButtonsAndText(self):
        self.T.delete('1.0', tk.END)
        selectMessage = """File Loaded.\nSelect Starting Molecule:"""
        self.T.insert(tk.END, selectMessage)
        self.addFileButton.pack_forget()
        self.singleCarbon = tk.Button(
            self.frame, text="Single Carbon", fg="black", bg='blue', command=self.buildModifyingState)
        self.singleCarbon.pack()
        self.randomMoleculeButton = tk.Button(
            self.frame, text="Random Molecule", fg="black", bg='blue', command=self.buildModifyingState)
        self.randomMoleculeButton.pack()

    def buildErrorState(self):
        self.T.delete('1.0', tk.END)
        errorMessage = """Incorrect File Selected.\nPlease Select Another File."""
        self.T.insert(tk.END, errorMessage)

    def buildModifyingState(self):
        self.env.reset()
        self.T.delete("1.0", tk.END)
        modifyingMessage = """Applying modifications to molecule...."""
        self.T.insert(tk.END, modifyingMessage)
        self.singleCarbon.pack_forget()
        self.randomMoleculeButton.pack_forget()
        # This function needs to loop until we determine the end...
        self.env.step()
        self.frame.after(2000)
        self.buildModifiedState()

    def buildModifiedState(self):
        modifiedMessage = """Modifications made to molecule."""
        self.T.delete("1.0", tk.END)
        self.T.insert(tk.END, modifiedMessage)
        self.renderButton = tk.Button(
            self.frame, text="Render Current Molecule", fg="black", bg='blue', command=self.renderMolecule)
        self.renderButton.pack()
        self.getModificationsListButton = tk.Button(
            self.frame, text="Get Modifications List", fg="black", bg='blue', command=self.produceModificationsList)
        self.getModificationsListButton.pack()

    def renderMolecule(self):
        self.env.render()
        self.renderButton.pack_forget()
        self.T.delete("1.0", tk.END)
        renderedMoleculeMessage = """Rendered Molecule"""
        self.T.insert(tk.END, renderedMoleculeMessage)
        self.label = Label(root, image=self.photo)
        self.label.pack()

    def restart_program(self):
        self.frame.destroy()
        self.frame = tk.Frame(root, bg="white")
        self.frame.pack(expand='true', fill="both")
        self.quitButton = tk.Button(self.frame, text="Quit", command=sys.exit)
        self.quitButton.pack(side='left')
        self.restartButton = tk.Button(
            self.frame, text="Restart", command=self.restart_program)
        self.restartButton.pack(side='right')
        self.T = Text(self.frame, height=2)
        Welcome = """Welcome! \nPlease Select a file to upload."""
        self.T.pack()
        self.T.insert(tk.END, Welcome)
        self.T.configure(font=("Courier", 16,))
        self.addFileButton = tk.Button(
            self.frame, text="Add File", fg="black", bg='blue', command=self.addFile)
        self.addFileButton.pack()
        self.label.destroy()

    def produceModificationsList(self):
        self.T.delete("1.0", tk.END)
        listPrinted = """Modifications List Produced as: output.txt"""
        self.T.insert(tk.END, listPrinted)

        os.system("open output.txt")


root = tk.Tk(className=' TEDD')
GuiWindow(root).pack()
root.mainloop()
