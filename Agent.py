#Importing Libraries
import random

#Class For Agent
class AgentClass:
    def agentAction(self,observationsList,reward):
        """
        Dummy Agent
        |||||||||||||||||||||||||||||||||||||||||||||||||||||||
        ||THIS CLASS IS FOR USER TO IMPLEMENT CODES FOR AGENT||
        |||||||||||||||||||||||||||||||||||||||||||||||||||||||
        Args:
            observationsList (list): List of observations returned from Step Function in Environment
            reward (float): Reward of Agent for its Previous Action

        Returns:
            networkx.Graph(): Returns Molecular Graph
        """        
        nextState = random.choice(observationsList)        
        return nextState