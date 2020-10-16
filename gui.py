#Importing Libraries
import gym
import random
import gym_molecule
import time
import sys
import multiprocessing
from rdkit import Chem
from kivy.app import App 
import kivy 
from tkinter import filedialog
from tkinter import Tk 
from kivy.uix.floatlayout import FloatLayout 
from tkinter.filedialog import askopenfilename
from kivy.uix.boxlayout import BoxLayout
from kivy.config import Config 
from kivy.uix.screenmanager import ScreenManager,Screen
from kivy.lang import Builder
from tkinter.messagebox import showinfo
from kivy.clock import Clock, mainthread
import atexit, os
from PIL import Image
from tkinter import messagebox
from pathlib import Path
import Agent

#Creating Config for Window Size
#This GUI works best with bigger Screens
Config.set('graphics', 'resizable', False)
Config.set('graphics', 'width', '1600')
Config.set('graphics', 'height', '800')

#Defining Global Variables
stop = True
startingMoleculeList = []
actionSpaceFileName = ""
actionSpaceFileUploaded = False
render = False
parallelThread = None
startingMoleculeFileUploaded = False
startingMoleculeFileName = ""
EnvLoops = None
StepsPerEnvLoop = None
renderer = None


class backgroundTab(Screen):

    def spinnerRenderingChange(self,text):
        """
        If Spinner is updated Global render value is also updated Accordingly
        Args:
            text (str): Text of what is Selected
        """        
        global render
        if text == "No Rendering":
            render = False
        else:
            render = True
    
    def spinnerAction(self,text):
        """
        Adds Button To GUI for User to upload a File
        Args:
            text (str): text of what is Selected
        """        
        if text == "Upload Action Space":
            self.buttonActionSpace.opacity=1
            self.buttonActionSpace.disabled= False
        else:
            self.buttonActionSpace.opacity=0
            self.buttonActionSpace.disabled= True
    
    def spinnerStartingMol(self,text):
        """
        Adds Button or InputText To GUI for User to upload a File or write in Starting Molecule

        Args:
            text (str): text of what is Selected
        """
        text = self.spinnerStartingMolecule.text 
        if text == "Random Molecule from List":
            self.buttonStartingMolecule.opacity=1
            self.buttonStartingMolecule.disabled= False
            self.inputForStartingMolecule.opacity=0
            self.inputForStartingMolecule.disabled= True
        elif text == "Certain Molecule":
            self.inputForStartingMolecule.opacity=1
            self.inputForStartingMolecule.disabled= False
            self.buttonStartingMolecule.opacity=0
            self.buttonStartingMolecule.disabled= True
        else:
            self.buttonStartingMolecule.opacity=0
            self.buttonStartingMolecule.disabled= True
            self.inputForStartingMolecule.opacity=0
            self.inputForStartingMolecule.disabled= True

    def actionSpaceFile(self):
        """
        Prompts User to upload a file containing the Action Space Chemical Reactions and handles errors and gives Warnings
        Returns:
            int 0 : Return
        """       
        global actionSpaceFileUploaded
        global actionSpaceFileName
        Tk().withdraw()
        filename = askopenfilename()
        my_file = Path(filename)
        if not my_file.is_file():
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Path Not Found For Action Space File")
            tempWindow.destroy()
            return 0
        if os.stat(filename).st_size == 0:
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Action Space File Is Empty")
            tempWindow.destroy()
            return 0

        actionSpaceFileUploaded = True
        actionSpaceFileName = filename
        return 0
            
        

    def startingMoleculeFile(self):
        """
        Prompts User to input a file containing list of Starting molecule for the starting state of the Environment and handles errors and gives Warnings
        Returns:
            int 0: return
        """     
        global startingMoleculeFileUploaded
        global startingMoleculeFileName
        Tk().withdraw()
        filename = askopenfilename()
        my_file = Path(filename)
        if not my_file.is_file():
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Path Not Found For Starting Molecule File")
            tempWindow.destroy()
            return 0
        if os.stat(filename).st_size == 0:
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", " Starting Molecule File Is Empty")
            tempWindow.destroy()
            return 0

        startingMoleculeFileUploaded = True
        startingMoleculeFileName = filename
        return 0
    
    def startEnvironment(self):
        """
        Method to control and determine the different variables to be used in the Environment such as the starting molecule, 
        the action space, Number of Environment Loops, Number of Steps per Environment Loop and Rendering
        Updates GUI by creating new widgets and removing old ones
        Creates and Starts running The Enviroment in parallel
        Returns:
            int 0: return
        """
        global actionSpaceFileName
        global render
        global actionSpaceFileUploaded
        global startingMoleculeList
        global EnvLoops
        global StepsPerEnvLoop

        #Gets the Starting Molecule or Molecules List, Handles different types of Errors
        if self.spinnerStartingMolecule.text == "Carbon":
            startingMoleculeList.append('C')
            
        elif self.spinnerStartingMolecule.text == "Random Molecule from List":
            if not startingMoleculeFileUploaded:
                tempWindow = Tk()
                tempWindow.withdraw()
                messagebox.showinfo("Warning", "Please upload Starting Molecules File")
                tempWindow.destroy()
                return 0
            else:
                file = open(startingMoleculeFileName)
                for lines in file.readlines():
                    try :
                        if Chem.MolFromSmiles(lines.strip("\n")):
                            startingMoleculeList.append(lines.strip("\n"))
                        else:
                            print("Invalid Starting Molecule")
                    except:
                        pass
                file.close()
        else:
            TempStartingMolecule = self.inputForStartingMolecule.text
            if TempStartingMolecule == "":
                tempWindow = Tk()
                tempWindow.withdraw()
                messagebox.showinfo("Warning", "Please Write in Starting Molecule")
                tempWindow.destroy()
                return 0
            elif not Chem.MolFromSmiles(TempStartingMolecule):
                tempWindow = Tk()
                tempWindow.withdraw()
                messagebox.showinfo("Warning", "Please Write Valid Starting Molecule")
                tempWindow.destroy()
                return 0
            startingMoleculeList.append(TempStartingMolecule)
            

        #Gets the Action Space, Handles different types of Errors
        if self.spinnerActionSpace.text == "Default":
            actionSpaceFileName = "DefaultActionSpace.txt"            
        else:
            if not actionSpaceFileUploaded:
                tempWindow = Tk()
                tempWindow.withdraw()
                messagebox.showinfo("Warning", "Please upload Action Space File")
                tempWindow.destroy()
                return 0

        #Gets the number of loops for environment and steps, Handles different types of Errors
        TempEnvLoops = self.inputEnvLoop.text
        TempStepsPerEnvLoop = self.inputEnvLoopStep.text
        if not TempEnvLoops.isdigit():
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Please Enter Positive Integer to number of Environment Loops")
            tempWindow.destroy()
            return 0
        else :
            EnvLoops = TempEnvLoops
        
        if not TempStepsPerEnvLoop.isdigit():
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Please Enter Positive Integer to Number of Steps per Environment Loop")
            tempWindow.destroy()
            return 0
        else :
            StepsPerEnvLoop = TempStepsPerEnvLoop

        #Checks if Environment Should be rendered, updates render accordingly
        if self.spinnerRendering.text == "No Rendering":
            render = False
        else:            
            render = True

        #Removes current widgets and adds new widgets
        self.buttonStart.opacity=0
        self.inputForStartingMolecule.opacity=0
        self.buttonStartingMolecule.opacity=0
        self.buttonActionSpace.opacity=0
        self.spinnerRendering.opacity=0
        self.spinnerStartingMolecule.opacity=0
        self.spinnerActionSpace.opacity=0
        self.lblStartingMolecule.opacity=0
        self.lblActionSpace.opacity=0
        self.buttonStart.disabled= True
        self.inputForStartingMolecule.disabled= True
        self.buttonStartingMolecule.disabled= True
        self.buttonActionSpace.disabled= True
        self.spinnerRendering.disabled= True
        self.spinnerStartingMolecule.disabled= True
        self.spinnerActionSpace.disabled= True
        self.lblStartingMolecule.disabled= True
        self.lblActionSpace.disabled= True
        self.lblRendering.opacity=0
        self.lblRendering.disabled= True
        self.outInfo.opacity=1
        self.outInfo.disabled= False
        self.outSettings.opacity=1
        self.outSettings.disabled= False
        if render:
            self.buttonRenderContinue.opacity = 1
            self.buttonRenderContinue.disabled = False
        self.buttonReset.opacity = 1
        self.buttonReset.disabled = False
        self.buttonHistory.opacity = 1
        self.buttonHistory.disabled = False
        self.imageMol.opacity = 1
        self.imageMol.disabled = False        
        self.buttonHistory.pos = (1360,90)
        self.inputEnvLoopStep.disabled= True 
        self.inputEnvLoop.disabled= True 
        self.inputEnvLoopStep.opacity= 0
        self.inputEnvLoop.opacity=0
        self.lblEnvironmentLoop.opacity=0
        self.lblPerEnvironmentLoop.opacity=0
        self.lblEnvironmentLoop.disabled= True 
        self.lblPerEnvironmentLoop.disabled= True 

        #For Optimization
        if not render:
            renderer.cancel()
            infoUpdate=Clock.schedule_interval(self.updateInfo, 0.5)

        #Starts the Environment in parallel to GUI thread
        self.environmentStepCaller()

        #Creates a Settings History Text
        tempText = ""
        if self.spinnerActionSpace.text == "Default":
            tempText = "Default Action Space was Chosen\n"
        else:
            tempText = "Action Space File was Uploaded.\n"
        
        if self.spinnerStartingMolecule.text == "Random Molecule from List":
            tempText += "Starting Molecule List was Uploaded\n"

        elif self.spinnerStartingMolecule.text == "Carbon":
            tempText += "Starting Molecule was Chosen as Carbon Atom\n"

        else:
            tempText += "Starting Molecule was Written in "+startingMoleculeList[0]+"\n"

        if render:
             tempText += "Rendering in 2D\n"
        else:
            tempText += "No Rendering was Chosen\n"
        tempText+=str(TempEnvLoops) + " Number of Environment Loops\n"
        tempText+=str(StepsPerEnvLoop) + " Number of Steps per Environment\n"
        self.outSettings.text = tempText       

        return 0

    
    def environmentStepCaller(self):
        """
        Creates a Thread for Gym Environment and starts it
        """        
        global parallelThread
        parallelThread=multiprocessing.Process(target = self.environmentStepCallerparallelThread,args=(render,startingMoleculeList,actionSpaceFileName,EnvLoops,StepsPerEnvLoop,))
        parallelThread.start()


    @staticmethod
    def environmentStepCallerparallelThread(render,startingMoleculeList,actionSpaceFileName,EnvLoops,StepsPerEnvLoop):
        """
        Thread which perform creates Enviroment and calls Steps from Environment and also 
        creates an instance of the Agent to learn and take actions on environment based on reward system

        Args:
            render (boolean): returns if render function should be enabled or disabled
            startingMoleculeList (list): list of molecules as the starting state of the environment
            actionSpaceFileName (string): Name of the Action Space file uploaded by the User or default
            EnvLoops (int): Number of Environment loops to be carried out as specified by the User
            StepsPerEnvLoop (int): Number of steps to be carried out on the states per Environment Cycle as specified by the User
        """
        #Starts The History of Enviroments
        history = "Starting The Environment\n"
        try:
            f=open("HistoryTemp.txt", "a")
            f.write(history)
            f.close()
        except:
            pass

        #Creates and initilizes The enviroment with help of Gym libraries
        env = gym.make('gym_molecule:molecule-v0')
        env.init(startingMoleculeList,actionSpaceFileName)

        #Creates instance of agent
        agent = Agent.AgentClass()

        #Starts Enviroment Loops, Reseting Everytime Step Function Returns Done as True
        for i_episode in range(int(EnvLoops)):   
            observation = []
            #Getting Starting Molecule by reseting Enviroment
            observation.append(env.reset())
            history="Starting Molecule " +str(Chem.MolToSmiles(env.currentState)) + "\n"
            history+="Environment loop number : "+str(i_episode+1)+"\n"
            try:
                f=open("HistoryTemp.txt", "a")
                f.write(history)
                f.close()
            except:
                pass
            #Taking First Step
            action = observation[0]
            observation, reward, done= env.step(action)

            #Starting Step loop
            for t in range(int(StepsPerEnvLoop)):
                if not len(observation)>0:
                    break
                action = agent.agentAction(observation,reward)
                observation, reward, done= env.step(action)
                #Adding History of Enviroment
                history=""
                if t>0:
                    history = "Step "+str(t)+"\n"
                    history+="Current State " +str(Chem.MolToSmiles(env.currentState)) + "\n"     
                history+="Reward For Previous Action " +str(reward) + "\n"
                
                #CAN BE CHANGED TO SPEED UP SYSTEM OR TO BE SLOWED DOWN FOR ANALYSIS
                time.sleep(1)

                if render:
                    env.render()

                #Checks if Step Function Returned Done as True and if it True Resets the Environment
                if done:
                    history+="Environment can't make anymore modifications to current state\n"
                    history+="Episode finished after {} timesteps".format(t+1) + "\n"
                    try:
                        f=open("HistoryTemp.txt", "a")
                        f.write(history)
                        f.close()
                    except:
                        pass
                    break

                try:
                    f=open("HistoryTemp.txt", "a")
                    f.write(history)
                    f.close()
                except:
                    pass
        #Closes The Enviroment
        env.close()


    def buttonRenderStopClicked(self):
        """
        Changse text on Button and enables function of Render button
        """        
        global render
        self.buttonRenderContinue.text = "Start Rendering"
        render = False

    def buttonRenderContinueClicked(self):
        """
        Changes text on Button and disable function of the Render Button
        """        
        global render
        self.buttonRenderContinue.text = "Stop Rendering"
        render = True

    def resetEnvironment(self):
        """
        Terminates the environment and gui and restarts afresh
        """        
        parallelThread.terminate()
        os.execl(sys.executable, os.path.abspath(__file__), *sys.argv)
    
    def exit(self):
        """
        Terminate the parallel threads
        Closes The GUI
        """        
        if not parallelThread == None:
            parallelThread.terminate()
        exit()
    
    def updateRender(self, *args):
        """
        Applies Rendering to the GUI by Updating Image Widget
        """        
        self.updateInfo()

        if render:
            try:
                self.imageMol.source = "RenderingImage.png"
                self.imageMol.reload()
            except:
                pass
        else:
            self.imageMol.source = "ImageFunny.jpg"
            self.imageMol.reload()
            pass

    def createHistory(self):
        """
        Deletes old History file and Makes New Empty one
        """        
        f=open("HistoryTemp.txt", "w")
        f.close()
        pass
    
    def updateInfo(self,*args):
        """
        Reads the contents in the History file and changes GUI accordingly
        """        
        try:
            f=open("HistoryTemp.txt", "r")
            contents = f.read()
            self.outInfo.text = contents
            f.close()
        except:
            pass
    
    def setRenderer(self,Renderer):
        """
        Sets global renderer by getting from argument

        Args:
            Renderer (kivy.Clock): Clock for updating GUI for Rendering
        """        
        global renderer
        renderer = Renderer

    def environmentHistory(self):
        """
        Saves the History file to the location selected by user
        Also Handles errors
        """        
        Tk().withdraw()
        filename = filedialog.askdirectory(initialdir='.')
        if os.path.isdir(filename):
            f=open("HistoryTemp.txt", "r")
            f2 = open(str(filename)+"\History.txt","w")
            text = f.read()
            f2.write(text)
            f2.close()
            f.close()
        else:
            tempWindow = Tk()
            tempWindow.withdraw()
            messagebox.showinfo("Warning", "Please Select Proper Directory")
            tempWindow.destroy()

    
class background(App):  
    def build(self):
        """
        Builds The GUI
        Returns:
            [background]: GUI Class
        """        
        im = Image.open("ImageFunny.jpg")
        im.save("RenderingImage.png")
        im.close()
        kv = Builder.load_file("background.kv")
        GUI = backgroundTab()
        GUI.createHistory()
        Renderer=Clock.schedule_interval(GUI.updateRender, 0.5)
        GUI.setRenderer(Renderer)
        return  GUI
    
    
def forceExit(text):
    """
    Terminates Threads that are running on background on Force Exit
    Args:
        text (str): exit
    """    
    if not parallelThread == None:
        parallelThread.terminate()


if __name__ == '__main__': 
    background().run()
    #On force Exit Calls forceExit Function to terminate Threads that are running on background
    atexit.register(forceExit,"exit")
    #For MultiThreading
    multiprocessing.freeze_support()