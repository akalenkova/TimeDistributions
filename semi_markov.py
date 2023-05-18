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
        self.transitions = transitions
        self.graph = {state : {"out_transitions" : {}, "in_transitions" : {}} for state in states}
        for transition in transitions :
            # print("Creating out_transition " + transition[0] + " -> " + transition[1])
            self.graph[transition[0]]["out_transitions"][transition[1]] = transition
            # print("Creating in_transition " + transition[1] + " -> " + transition[0])
            self.graph[transition[1]]["in_transitions"][transition[0]] = transition
    
    def get_cheapest_node(self,states) : 
        cheapest_state = ""
        cheapest_cost = float("inf")
        for state in states :
            cost = self.get_computation_cost(state)
            if cost < cheapest_cost :
                cheapest_cost = cost
                cheapest_state = state
        return cheapest_state

    def reduce_all(self, states, log):
        state_size = len(states)
        for i in range(1,state_size) :
            cheapest_node = self.get_cheapest_node(states)
            states.remove(cheapest_node)
            label = "State " + str(i) + " out of " + str(state_size)
            self.reduce_node(cheapest_node, label, log)
            # if(cheapest_node != 'start' and cheapest_node != 'end'):
                # print("Deleting cheapest node : " + cheapest_node)
                # self.graph.pop(cheapest_node)     

    def reduce_node(self, state, label, log):
        if ((state == 'start') or (state == 'end')):
            return
        
        # print("\n")
        # print("Removing node " + state)
        # Calsulate self-loop time
        self_loop_time = MultiGauss([1], [Gauss(0,0)])
        self_loops = set()

        self_loop_time = self.calculate_self_loop_time(state, 0.0001)

        #  Add new transitions
        in_transitions = self.get_in_transitions(state)
        out_transitions = self.get_out_transitions(state)

        # print("Number of in_transitions of node " + state + " : " + str(len(in_transitions)))
        # print("Number of out_transition of node " + state + " : " + str(len(out_transitions)))
        count = 0

        # self.print_graph()

        for in_state, v in list(in_transitions.items()):
            for out_state, vv in list(out_transitions.items()):
                # print("In state : " + in_state)
                # print("Out state : " + out_state)
                print(count)
                count+= 1
                if in_state != out_state:
                    # print("in_state : " + str(v))
                    # print("out_state : " + str(vv))

                    p = self.get_probability(in_state, out_state)
                    time = self.get_time(in_state, out_state)

                    in_state_to_state_prob = self.get_probability(in_state,state)
                    state_to_out_state_prob = self.get_probability(state,out_state)
                    state_to_state_prob = self.get_probability(state,state)

                    # print("in_state_to_state_prob : " + str(in_state_to_state_prob))
                    # print("state_to_out_state_prob : " + str(state_to_out_state_prob))
                    # print("state_to_state_prob : " + str(state_to_state_prob))

                    new_prob = in_state_to_state_prob * state_to_out_state_prob/(1-state_to_state_prob)

                    all_p = p + new_prob

                    # print("p : " + str(p))
                    # print("new_prob : " + str(new_prob))
                    # print("all_p : " + str(all_p))
                    m1 = self.get_time(in_state, state)
                    m2 = self.get_time(state, out_state)
                    new_time = mult_gauss_convolution(m1,self_loop_time)
                    new_time = mult_gauss_convolution(new_time, m2)
                    all_time = mult_gauss_sum(time, new_time, p/all_p, new_prob/all_p)

                    new_transition = (in_state, out_state, all_p, all_time)
                    # Add new transition
                    self.add_out_state_transition(in_state,out_state,new_transition)
                    self.add_in_state_transition(out_state,in_state,new_transition)
                    mean_log_time = mean_time_between_events(in_state,out_state,[state],log)
                    # print("Removing out_state")
                    self.remove_transition(state, out_state)
            # print("Removing in_state")
            # removing the inwards and outwards transitions
            self.remove_transition(in_state, state)

    # def print_graph(self):
        # for state, item in self.graph.items() :
            # print("State : " + state)
            # print("Out Transitions : ", end="")
            # for out_transition in self.graph[state]["out_transitions"]:
                # print(" " + out_transition + "|", end="")
            # print("")
            # print("In Transitions : ", end="")
            # for out_transition in self.graph[state]["in_transitions"]:
                # print(" " + out_transition + "|", end="")
            # print("\n")


    def remove_transition(self, start, end):
        # print("Trying to remove " +  start + "->" + end)
        if start != "start" and start != "end" and end != "start" and end != "end" :
            self.graph[start]["out_transitions"].pop(end)
            # print("Removing out transition from " + start + " -> " +  end)
        if start != "start" and start != "end" and end != "start" and end != "end" :
            self.graph[end]["in_transitions"].pop(start)
            # print("Removing in transition from " + end + " -> " +  start)

    def add_out_state_transition(self, from_state, to_state, new_transition):
        # print("Adding new out_state from " + from_state + " -> " + to_state)
        self.graph[from_state]["out_transitions"][to_state] = new_transition

    def add_in_state_transition(self, from_state, to_state ,new_transition):
        # print("Adding new in_state from " + from_state + " -> " + to_state)
        self.graph[from_state]["in_transitions"][to_state] = new_transition

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