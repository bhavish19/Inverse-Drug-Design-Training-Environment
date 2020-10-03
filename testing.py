import gym 
import gym_molecule

#* Used for testing our molecule enviroment with corresponding functions
#* This will be created in the GUI
env = gym.make("gym_molecule:molecule-v0")

env.reset()

env.step()