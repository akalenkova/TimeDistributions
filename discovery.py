from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.objects.dfg.exporter import exporter as dfg_exporter
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery
from pm4py.statistics.start_activities.log import get as start_activities
from pm4py.statistics.end_activities.log import get as end_activities
from pm4py.visualization.dfg import visualizer as dfg_visualization
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import statsmodels.api as sm
from itertools import groupby #cant find
import gamma_distribution
import numpy as np
import seaborn as sns
from fit_distribution import fit_gauss, find_peaks_lib
import collections
import log_parser
from pm4py import view_dfg, discover_dfg
import time 
from copy import deepcopy
from mult_gauss import MultiGauss
from gauss import Gauss
from semi_markov import SemiMarkov
import heapq


def extract_times():
    for trace in log:
        first = True
        for next_event in trace:
            if not first:
                time = next_event['time:timestamp'] - event['time:timestamp']
                hours = time.total_seconds()//3600
                if hours > 250:
                    continue
                if not event['concept:name'] in times_dictionary.keys():
                    times_dictionary[event['concept:name']] = [time.total_seconds()//3600]
                else:
                    times_dictionary[event['concept:name']].append(time.total_seconds()//3600)
            event = next_event
            first = False

def extract_times_with_future(log):
    for trace in log:
        first = True
        for next_event in trace:
            if not first and next_event['concept:name'] != 'end' and event['concept:name'] != 'start':
                time = next_event['time:timestamp'] - event['time:timestamp']
                if not event['concept:name'] + '->' + next_event['concept:name'] in times_dictionary.keys():
                    times_dictionary[event['concept:name'] + '->' + next_event['concept:name']] = [time.total_seconds()//3600]
                else:
                    times_dictionary[event['concept:name'] + '->' + next_event['concept:name']].append(time.total_seconds()//3600)
            event = next_event
            first = False

def extract_times_event_log():
    times = []
    for trace in log:
        start = trace[0]['time:timestamp']
        end = trace[len(trace)-1] ['time:timestamp']   
        time = end - start
        times.append(time.total_seconds()//3600)
    #print(times)
    return times

# retrieve distribution of values from given y
def retrieve_distribution(y):
    result = {}
    for i in y:
        count_i = y.count(i)
        result[i] = count_i
        result[i] /= len(y)
    #gamma_distribution.fit_gamma(result)
    return result

def build_semi_markov(dfg, multi_gausses):
    
    states = set()
    transitions = set()
    out_frequences = {}

    for key in dfg.keys():
        states.add(key[0])
        #print(states)
    
    for key1 in states:
        out_frequences[key1] = 0
        for key2 in states:
            if dfg[key1,key2] > 0:
                out_frequences[key1] += dfg[key1,key2]


    for key1 in states:
        for key2 in states:
            if dfg[key1,key2] > 0:
                if ((key1 == 'start') or  (key1 == 'end') or (key2 == 'start') or  (key2 == 'end')):
                    transitions.add(tuple([key1, key2, dfg[key1,key2]/out_frequences[key1], MultiGauss([1], [Gauss(0, 0)])]))
                else:
                    transitions.add(tuple([key1, key2, dfg[key1,key2]/out_frequences[key1], 
                    multi_gausses["['" + str(key1) + "', '" + str(key2) + "']"]]))
    print('DFG is built')
    time.sleep(10)
    return SemiMarkov(states, transitions)


variant = xes_importer.Variants.ITERPARSE
parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
log = xes_importer.apply('/Users/alexhu/Documents/Github/TimeDistributions/logs/bpi_challenge_2013_incidents.xes', 
    variant=variant, parameters=parameters)
for k in [1,2,3,4]:
    start = time.time()
    log_for_discovery = deepcopy(log)
    times_dictionary = {}
    log_processed = log_parser.prepare_log(log_for_discovery, k)

    # the initial discovery
    dfg, start_activities, end_activities = discover_dfg(log_processed)
        #for s in start_activities.keys():
        #    start_activities[s] = 0
    dfg["end", "start"] = 1
    #view_dfg(dfg, start_activities, end_activities)

    extract_times_with_future(log_processed)


    """
    Fitting using Gaussian KDE
    """
    # enchancement
    mult_gausses = {}
    #fig = plt.figure(figsize=(10,220))
    i = 1
    for key in sorted(times_dictionary.keys()):
    #  if i<9:
        i += 1
    #      continue
        #print('Activity:')
        #print(key)

        #dist = retrieve_distribution(times_dictionary.get(key))
        #h = fig.add_subplot(len(times_dictionary), 1, i)
        #y, x, _ = plt.hist(times_dictionary.get(key), bins=200, fc=(0, 1, 0, 0.4), density = True, label='Histogram for time distribution')
        #print(x)
        #print(y)
        #plt.plot(x, y, 'tab:brown', linewidth=1)
        #h.set_ylim(0, min(max(y),1))
        #h.title.set_text(key)
        #h.set_xlabel('Waiting time in hours')
        #h.set_ylabel("Number of occurancies")
        #print(times_dictionary.get(key))
        kde = sm.nonparametric.KDEUnivariate(times_dictionary.get(key))
        kde.fit(bw=4, kernel='gau')  # Estimate the densities
    

        #print(kde.cdf)
        #print(len(kde.density))
        #print(len(kde.support))
        #h.plot(kde.support, kde.density, label="KDE")
        multi_gauss = fit_gauss(kde.support, kde.density, key)
        mult_gausses[str([key.partition('->')[0], key.partition('->')[2]])] = multi_gauss

        #y = retrieve_distribution(times_dictionary.get(key))
        #print(y) 
        #od = collections.OrderedDict(sorted(y.items()))
        #fit_gauss(list(od.keys()), list(od.values()))
        
        #print(od.values())
        #fit_gauss(od.keys(), od.values())

        #print(mult_gausses)
        i += 1

    # builds the new markov chain
    semi_markov = build_semi_markov(dfg, mult_gausses)
    states = deepcopy(semi_markov.states)
    i = 1


    semi_markov.reduce_all(states, log)

    # for state in states:
    #     #if str(state) != "Queued":
    #     #    continue
    #     # print("states")
    #     # print(states)
    #     # print("transitions")
    #     # print(semi_markov.transitions)
    #     # print("\n\n\n")
    #     label = "State " + str(i) + " out of " + str(len(states))

    #     semi_markov.reduce_node(state, label, log)
    #     i += 1  

    # outputs the results
    print("Number of transitions: " + str(len(semi_markov.transitions)))
    for transition in semi_markov.transitions:
        if transition[0] == 'start':
            print(transition)
            multi_gauss = transition[3]
            multi_gauss.remove_zero()
            color = {
                1: "tab:red",
                2: "tab:blue",
                3: "tab:orange",
                4: "tab:brown",
                5: "tab:purple",
            }
            multi_gauss.plot_mult_gauss(range(-10,1000,1), label="Semi-Markov Model k="+str(k), color = color.get(k))
            print()
            print('Overall: ---------------')
            print(multi_gauss.calculate_mean())

    end = time.time()
    print("Execution time for k=" + str(k))
    print(end - start)

"""
Approxiamting event log
"""
event_log_times = extract_times_event_log()
filtered_event_log_times = []
for time in event_log_times:
    if time < 1000:
        filtered_event_log_times.append(time)
#print(event_log_times)
y, x, _ = plt.hist(filtered_event_log_times, bins=100, fc=(0, 1, 0, 0.4), density=True, label='Event log')
#plt.xlabel('Waiting time in hours')
#plt.ylabel('Probability')
#kde_log = sm.nonparametric.KDEUnivariate(list(event_log_times))
#kde_log.fit(bw=20, kernel='gau')  # Estimate the densities
#plt.title('Event log time')
#multi_gauss_log = fit_gauss(kde_log.support, kde_log.density, 'Event log')
plt.xlim([-10, 1000])
plt.legend(loc="upper right")
plt.title('')
plt.show()

#plt.savefig('/Users/a1230101//Documents/GitHub/TimeDistributions/time_plots/DomesticDeclarations.pdf')

"""
Fitting without KDE
"""
"""
fig = plt.figure(figsize=(10,220))
i = 1
for key in sorted(times_dictionary.keys()):
    times = times_dictionary.get(key)
    h = fig.add_subplot(len(times_dictionary), 1, i)
    y_hist, x_hist, _ = h.hist(times, bins=1000, facecolor='g', density = True,)
    h.set_ylim(0, min(max(y_hist),1))
    h.title.set_text(key)
    h.set_xlabel('Waiting time in hours')
    h.set_ylabel("Number of occurancies")
    dict_of_times = retrieve_distribution(times)
    odered_times = collections.OrderedDict(sorted(dict_of_times.items()))
    fit_gauss(x_hist[1:], y_hist)

    #h.plot(kde.support, kde.density, label="KDE")
      
    i += 1
#plt.savefig('/Users/a1230101//Documents/GitHub/TimeDistributions/time_plots/DomesticDeclarationsFit.pdf')
"""