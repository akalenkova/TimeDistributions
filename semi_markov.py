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

class Graph: 
        

class SemiMarkov:
    def __init__(self, states, transitions):
        self.graph = {state : {"out_transitions" : {}, "in_transitions" : {}} for state in states}
        for transition in transitions :
            print("Creating edge " + transition[0] + " -> " + transition[1])
            self.adjacency_list[transition[0]]["out_transitions"][transition[1]] = transition
            self.adjacency_list[transition[1]]["in_transitions"][transition[0]] = transition
    
    def get_cheapest_node(self) : 
        cheapest_state = ""
        cheapest_cost = float('inf')
        for state in states :
            cost = self.get_computation_cost(state)
            if cost < cheapest_cost :
                cheapest_cost = cost
                cheapest_state = 
                
        return cheapest_state

    def reduce_all(self, states, log):

        for i in range(1,len(states)) :
            cheapest_node = get_cheapest_node()
            label = "State " + str(i) + " out of " + str(len(states))
            self.reduce_node(cheapest_node, label, log)
        

    def reduce_node(self, state, label, log):
        if ((state == 'start') or (state == 'end')):
            return
        else:
            print("Removing node " + state)
            # Calsulate self-loop time
            self_loop_time = MultiGauss([1], [Gauss(0,0)])
            self_loops = set()

            self_loop_time = self.calculate_self_loop_time(state, 0.0001)

            #  Add new transitions
            in_transitions = self.get_in_transitions(state)
            out_transitions = self.get_out_transitions(state)

            for in_state, v in list(in_transitions.items()):
                for out_state, vv in list(out_transitions.items()):
                    if in_state != out_state:
                        # print(in_state)
                        # print(out_state)
                        # print("in_state : " + str(v))
                        # print("out_state : " + str(vv))

                        p = self.get_probability(in_state, out_state)
                        time = self.get_time(in_state, out_state)
                        new_p = self.get_probability(in_state,state)*self.get_probability(state,out_state)/(1-self.get_probability(state,state))
                        all_p = p + new_p

                        print("p : " + str(p))
                        print("new_p : " + str(new_p))
                        print("all_p : " + str(all_p))
                        m1 = self.get_time(in_state, state)
                        m2 = self.get_time(state, out_state)
                        new_time = mult_gauss_convolution(m1,self_loop_time)
                        new_time = mult_gauss_convolution(new_time, m2)
                        all_time = mult_gauss_sum(time, new_time, p/all_p, new_p/all_p)

                        # removing the inwards and outwards transitions
                        self.remove_out_state_node(out_state, state)
                        self.remove_in_state_node(in_state, state)

                        # Add new transition
                        transition = (in_state, out_state, all_p, all_time)
                        self.add_out_state_transition(in_state,out_state,transition)
                        self.add_in_state_transition(out_state,in_state,transition)
                        mean_log_time = mean_time_between_events(in_state,out_state,[state],log)
            self.graph.pop(state)

    def remove_out_state_node(self, to_node, remove_node):
        self.graph[to_node]["in_transition"].pop(remove_node)

    def remove_in_state_node(self, from_node, remove_node):
        self.graph[from_node]["out_transitions"].pop(remove_node)

    def add_out_state_transition(self, from_state, to_state, transition):
        if(to_state in self.graph[out_state]["out_transitions"] ):
            # perform merge
        else :
            self.graph[from_state]["out_transitions"][to_state] = transition

    def add_in_state_transition(self, from_state, to_state , transition):
        if(to_state in self.graph[out_state]["in_transition"]):
            #perform merge
        self.graph[from_state]["in_transition"][to_state] = transition

    def  calculate_self_loop_time(self, state, threshold):
        m1 = self.get_time(state, state)
        p = self.get_probability(state, state)
        m = MultiGauss([1-p],[Gauss(0,0)])
        p_current = p * (1-p)
        k = 1
        conv = MultiGauss([1],[Gauss(0,0)])

        while (p_current > threshold):
            conv = mult_gauss_convolution(m1, conv)
            m = mult_gauss_sum(m, conv, 1, p_current)
            p_current *= p
            k += 1
        return m
    
    def get_in_transitions(self, state):
        # for key, value in self.graph[state]["in_transitions"].items() :
        #     print (key, value)
        return self.graph[state]["in_transitions"]

    # change
    def get_self_loop(self, state):
        return self.graph[state]["in_transitions"]
    
    def get_out_transitions(self, state):
        # for key, value in self.graph[state]["out_transitions"].items() :
        #     print (key, value)
        return self.graph[state]["out_transitions"]
        
    def get_computation_cost(self, state):
        return len(self.get_in_transitions(state)) * len(self.get_out_transitions(state)) 
    
 
    def get_probability(self, state1, state2):
      
        if (state1 in self.graph) and (state2 in self.graph[state1]["out_transitions"]):
            return self.graph[state1]["out_transitions"][state2][2]
        return 0


    def get_time(self, state1, state2):
        if (state1 in self.graph) and (state2 in self.graph[state1]["out_transitions"]):
            return self.graph[state1]["out_transitions"][state2][3]
        return MultiGauss([],[])