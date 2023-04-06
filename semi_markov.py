import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special
from mult_gauss import MultiGauss
from gauss import Gauss
import math
from convolution import mult_gauss_convolution, mult_gauss_sum, mult_gauss_self_convolution

def mean_time_between_events(e1,e2,skip_events,log):
        times = set()
        for trace in log:
            for i in range(len(trace) - 1):
                if trace[i]['concept:name'] == e1:
                    for k in range (i+1, len(trace)):
                        if trace[k]['concept:name'] == e2:
                            time = trace[k]['time:timestamp'] - trace[i]['time:timestamp']
                            times.add(time.total_seconds()//3600)
                            break
                        if not (trace[k]['concept:name'] in skip_events):
                            break
        if len(times) > 0:
            mean = sum(times)/len(times)
        else: 
            mean = 0
        return mean

class SemiMarkov:
    def __init__(self, states, transitions):
        self.states = states 
        # Each transition is a tuple: (from, to, probability, multi Gauss)
        self.transitions = transitions 

    

    # This is the code that I need to touch 
    def reduce_node(self, state, label, log):
        #print("Removing state " + str(state))
        if ((state == 'start') or (state == 'end')):
            return
        else:
            # Calsulate self-loop time
            self_loop_time = MultiGauss([1], [Gauss(0,0)])
            #print("Start calculating self-loop time...")
            self_loops = set()
            for transition in self.transitions:
                if ((transition[0] == state) and (transition[1] == state)):
                    self_loops.add(transition)
            for transition in self_loops:
                self_loop_time = self.calculate_self_loop_time(state, 0.0001)
#                self_loop_time.plot_mult_gauss(range(0,1000,1), str(state) + " " + str(state))
                #print("Self-loop mean:")
                #print(self_loop_time.calculate_mean())
                #plt.title(label=str(state) + " " + str(state))
                #plt.show()
            #print("End calculating self-loop time...")
            #  Add new transitions
            in_transitions = self.get_in_transitions(state)
            out_transitions = self.get_out_transitions(state)
            number_of_new_transitions = len(in_transitions) * len(out_transitions)
            i = 1
            for in_transition in in_transitions:
                in_state = in_transition[0]
                for out_transition in out_transitions:
                    out_state = out_transition[1]
                    if in_state != out_state:
                        #print("Adding transiton from " + str(in_state) + " to " + str(out_state))
                        i += 1
                        p = self.get_probability(in_state, out_state)
                        time = self.get_time(in_state, out_state)
                        #print("In probability:")
                        #print(self.get_probability(in_state,state))
                        #print("Out probability:")
                        #print(self.get_probability(state,out_state))
                        #print("Self probability:")
                        #print(self.get_probability(state,state))
                        new_p = self.get_probability(in_state,state)*self.get_probability(state,out_state)/(1-self.get_probability(state,state))
                        all_p = p + new_p
                        m1 = self.get_time(in_state, state)
                        #print("First mean time:")
                        #print(m1.calculate_mean())
                        m2 = self.get_time(state, out_state)
                        #print("Second mean time:")
                        #print(m2.calculate_mean())
                        new_time = mult_gauss_convolution(m1,self_loop_time)
                        new_time = mult_gauss_convolution(new_time, m2)
                        #print("Old time:")
                        #print(time.calculate_mean())
                        #print("New time:")
                        #print(new_time.calculate_mean())
                        #print("Old probaility:")
                        #print(p)
                        #print("New probaility:")
                        #print(new_p)
                        all_time = mult_gauss_sum(time, new_time, p/all_p, new_p/all_p)
                        #print("All mean time:")
                        #print(all_time.calculate_mean())

                        # Remove old transition
                        transition_to_remove = set()
                        for transition in self.transitions:
                            if ((transition[0] == in_state) and (transition[1] == out_state)):
                                transition_to_remove.add(transition)
                        for transition in transition_to_remove:
                            self.transitions.remove(transition)
                            
                        # Add new transition
                        self.transitions.add((in_state, out_state, all_p, all_time))

#                        all_time.plot_mult_gauss(range(0,1000,1), str(in_state) + " " + str(out_state))
                        #print("Sum of probabilitties:")
                        #print(all_time.calculate_sum_probabilities())
                        #print("Mean:")
                        #print(all_time.calculate_mean())
                        mean_log_time = mean_time_between_events(in_state,out_state,[state],log)
                        #print(mean_log_time)
                        #plt.title(label=str(in_state) + " " + str(out_state))
                        #plt.show()

            # Remove state
            transition_to_remove = set()
            for transition in self.transitions:
                if ((transition[0] == state) or (transition[1] == state)):
                    transition_to_remove.add(transition)
            for transition in transition_to_remove:
                 self.transitions.remove(transition)
            self.states.remove(state)


    def  calculate_self_loop_time(self, state, threshold):
        m1 = self.get_time(state, state)
        p = self.get_probability(state, state)
        #print("Probability of the self-loop " + str(p))
        m = MultiGauss([1-p],[Gauss(0,0)])
        p_current = p * (1-p)
        k = 1
        conv = MultiGauss([1],[Gauss(0,0)])
        while (p_current > threshold):
            #print(p_current)
            #print("Calculate self-loop:")
            #print(k)
            conv = mult_gauss_convolution(m1, conv)
            m = mult_gauss_sum(m, conv, 1, p_current)
            #print(m.calculate_mean())
            p_current *= p
            k += 1
        #print(m.probabilities)
        #print(m.gaussians)
        #print("---------------------------------")
        return m
    
    def get_in_transitions(self, state):
        in_transitions = set()
        for transition in self.transitions:
            if (transition[1] == state) and (transition[0] != state):
                in_transitions.add(transition)
        return in_transitions
    
    def get_out_transitions(self, state):
        out_transitions = set()
        for transition in self.transitions:
            if (transition[0] == state) and (transition[1] != state):
                out_transitions.add(transition)
        return out_transitions
 
    def get_probability(self, state1, state2):
        for transition in self.transitions:
            if ((transition[0] == state1) and (transition[1] == state2)):
                return transition[2]
        return 0
    
    def get_time(self, state1, state2):
        for transition in self.transitions:
            if ((transition[0] == state1) and (transition[1] == state2)):
                return transition[3]
        return MultiGauss([],[])
    
    


    