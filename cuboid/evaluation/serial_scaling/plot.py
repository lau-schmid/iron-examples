#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
import copy
from sets import Set

# load format package
SCENARIO='cuboid'
#sys.path.insert(0, src_path+"/evaluation")
import format as fo

# determine if plots are shown
show_plots = True
if len(sys.argv) >= 2:
  show_plots = False
  
remove_outlier = False
outlier_top = 1
outlier_bottom = 1
  
# read csv file
report_filename = "duration_timing.csv";
print "report file: {}".format(report_filename)
data = []
with open(report_filename) as csvfile:
  reader = csv.reader(csvfile, delimiter=';')
  for row in reader:
    if len(row) > 0:
      if '#' not in row[0]:
        data.append(row)

n = len(data)

def getCol(colno):
  cols = []
  for i in range(len(data)):
    cols.append(data[i][colno])
  return cols  
  
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False
    
# 0  Stamp
# 1  Host
# 2  NProc
# 3  X
# 4  Y
# 5  Z
# 6  F
# 7  Total FE
# 8  Total M
# 9  End Time
# 10 Dur. Init
# 11 Stretch Sim
# 12 Int. Init
# 13 Main Sim
# 14 Total
# 15 Total (User)
# 16 Total (System)
# 17 ODE
# 18 Parabolic
# 19 FE
# 20 FE before Main Sim

max_index = 20
int_indices = [2, 3, 4, 5, 6, 7, 8, 9]
float_indices = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

def extract_data(data):

  datasets = dict()
  
  # extract lines
  n = 0
  for dataset in data:
    
    if len(dataset) < 17:
      print "Warning: invalid data set"
      continue
    
    # copy dataset to new_data
    new_data = dataset
    
    # extract some values
    
    for index in int_indices:
      new_data[index] = int(new_data[index])     if isint(new_data[index])    else 0
    for index in float_indices:
      new_data[index] = float(new_data[index])     if isfloat(new_data[index])    else 0.0
      
    # define sorting key
    key = "{:03d}".format(new_data[7])+"|"+str(new_data[6])  # sort by FE
    #key = "{:06d}".format(new_data[8])  # sort by M
      
    # fill dataset to size of max_index
    if len(dataset) < max_index:
      a += (len(dataset)-max_index) * [0]
      
    # store extracted values
    if key not in datasets:
      datasets[key] = dict()
      datasets[key]["value"] = []
    
    datasets[key]['value'].append(new_data)
    n += 1
    
  # compute mean value
  for key in datasets:
    
    result = copy.copy(datasets[key]['value'][0])
    variances = copy.copy(datasets[key]['value'][0])
    result += [0]
    
    for i in int_indices + float_indices:
      
      # reduce values
      result[i] = 0
      value_list = []
      for j in range(len(datasets[key]['value'])):
        value = datasets[key]['value'][j][i]
        
        if value != 0:
          value_list.append(value)
      
      # remove outlier
      value_list = sorted(value_list)
      n = len(value_list)
      
      #print "i={}, value_list for {}: {}".format(i, key, value_list)
      
      if n > outlier_bottom+outlier_top and remove_outlier:
        value_list = value_list[outlier_bottom:-outlier_top]
        
      # compute mean and standard deviation
      result[i] = np.mean(value_list)
      variances[i] = np.std(value_list)
      #print "mean: {}, std: {}".format(result[i], variances[i])
        
    result[max_index+1] = len(datasets[key]["value"])
        
    datasets[key]['value'] = result
    datasets[key]['variance'] = variances
    
  datasets = collections.OrderedDict(sorted(datasets.items()))
    
  return datasets
  
datasets = extract_data(data)
#######################################################
# plot
# x-axis: fe
# y-axis: user time
# datasets: f

plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 2
output_path = ""
plotdata = dict()
xdata = Set()
plotkeys = Set()

for key in datasets:
  
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  f = dataset[6]
  fe = dataset[7]
  user = dataset[15]
  ode_time = dataset[17]
  parabolic_time = dataset[18]
  fe_time = dataset[19]
  pre_fe_time = dataset[20]
  
  plotkey = (f)
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = dict()
    plotdata[plotkey]['variance'] = dict()
    
  xvalue = fe
  yvalue = user
  yvalue_variance = 0
    
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)

for plotkey in plotkeys:
  print plotkey
  
xlist = list(xdata)

######################
# plot times
plt.figure(2, figsize=(10,8))
plotdata = dict()
xdata = Set()
plotkeys = Set()

colors = {
  (1) : "ro",
  (2) : "go",
  (3) : "bo",
  (4) : "co",
  (5) : "mo",
}
labels = {
  (1) : "F=1",
  (2) : "F=2",
  (3) : "F=3",
  (4) : "F=4",
  (5) : "F=5",
}
for key in datasets:
  
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  f = dataset[6]
  fe = dataset[7]
  user = dataset[14]
  
  ode_time = dataset[17]
  parabolic_time = dataset[18]
  fe_time = dataset[19]
  pre_fe_time = dataset[20]
  
  ode_time_v = variances[17]
  parabolic_time_v = variances[18]
  fe_time_v = variances[19]
  pre_fe_time_v = variances[20]
  
  #print "f=",f,", fe=",fe
  if f != 1.0:
    continue
  
  ## ode
  plotkey = 0
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = dict()
    plotdata[plotkey]['variance'] = dict()
    
  xvalue = fe
  yvalue = ode_time
  yvalue_variance = ode_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## parabolic
  plotkey = 1
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = dict()
    plotdata[plotkey]['variance'] = dict()
    
  xvalue = fe
  yvalue = parabolic_time
  yvalue_variance = parabolic_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## fe_time
  plotkey = 2
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = dict()
    plotdata[plotkey]['variance'] = dict()
    
  xvalue = fe
  yvalue = fe_time
  yvalue_variance = fe_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## pre_fe
  plotkey = 3
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = dict()
    plotdata[plotkey]['variance'] = dict()
    
  xvalue = fe
  yvalue = pre_fe_time
  yvalue_variance = pre_fe_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
xlist = list(xdata)

colors = {
  0 : "yo-",
  1 : "ro-",
  2 : "go-",
  3 : "bo-",
}
labels = {
  0 : "Duration ODE solver",
  1 : "Duration Parabolic solver",
  2 : "Duration FE solver (main sim)",
  3 : "Duration FE solver (pre-sim)",
}

for plotkey in plotkeys:
  xlist = list(plotdata[plotkey]["value"])
  ylist = [y for y in plotdata[plotkey]["value"].values()]
  yerr = [y for y in plotdata[plotkey]['variance'].values()]

  plt.errorbar(xlist, ylist, fmt=colors[plotkey], yerr=yerr, label=labels[plotkey])
  
ax = plt.gca()
#ax.set_yscale('log', basey=10) 
plt.title(u'duration solvers, cuboid muscle, serial run, F=1')
plt.xlabel('finite elements total')
plt.ylabel('duration (s)')
plt.legend(loc='best')
plt.grid(which='both')
plt.tight_layout()
plt.savefig(output_path+SCENARIO+'_timing_f1.png')

######################
# plot times F=2
plt.figure(3, figsize=(10,8))
plotdata = dict()
xdata = Set()
plotkeys = Set()

for key in datasets:
  
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  f = dataset[6]
  fe = dataset[7]
  user = dataset[14]
  
  ode_time = dataset[17]
  parabolic_time = dataset[18]
  fe_time = dataset[19]
  pre_fe_time = dataset[20]
  
  ode_time_v = variances[17]
  parabolic_time_v = variances[18]
  fe_time_v = variances[19]
  pre_fe_time_v = variances[20]
  
  #print "f=",f,", fe=",fe
  if f != 2.0:
    continue
  
  ## ode
  plotkey = 0
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  xvalue = fe
  yvalue = ode_time
  yvalue_variance = ode_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## parabolic
  plotkey = 1
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  xvalue = fe
  yvalue = parabolic_time
  yvalue_variance = parabolic_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## fe_time
  plotkey = 2
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  xvalue = fe
  yvalue = fe_time
  yvalue_variance = fe_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## pre_fe
  plotkey = 3
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  xvalue = fe
  yvalue = pre_fe_time
  yvalue_variance = pre_fe_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
xlist = list(xdata)

for plotkey in plotkeys:
  xlist = list(plotdata[plotkey]["value"])
  ylist = [y for y in plotdata[plotkey]["value"].values()]
  yerr = [y for y in plotdata[plotkey]['variance'].values()]

  plt.errorbar(xlist, ylist, fmt=colors[plotkey], yerr=yerr, label=labels[plotkey])
  
ax = plt.gca()
#ax.set_yscale('log', basey=10) 
plt.title(u'Laufzeit F=2')
plt.xlabel('Anzahl Finite Elemente')
plt.ylabel('Laufzeit (user) (s)')
plt.legend(loc='best')
plt.grid(which='both')
plt.tight_layout()
plt.savefig(output_path+SCENARIO+'_timing_f2.png')

######################
# plot times 1D elements
plt.figure(4, figsize=(10,8))
plotdata = dict()
xdata = Set()
plotkeys = Set()

for key in datasets:
  
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  f = dataset[6]
  fe = dataset[7]
  user = dataset[14]
  m = dataset[8]
  
  ode_time = dataset[17]
  parabolic_time = dataset[18]
  fe_time = dataset[19]
  pre_fe_time = dataset[20]
  
  ode_time_v = variances[17]
  parabolic_time_v = variances[18]
  fe_time_v = variances[19]
  pre_fe_time_v = variances[20]
  
  ## ode
  plotkey = 0
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  xvalue = m
  yvalue = ode_time
  yvalue_variance = ode_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## parabolic
  plotkey = 1
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  yvalue = parabolic_time
  yvalue_variance = parabolic_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## fe_time
  plotkey = 2
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  yvalue = fe_time
  yvalue_variance = fe_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
  ## pre_fe
  plotkey = 3
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  yvalue = pre_fe_time
  yvalue_variance = pre_fe_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
xlist = list(xdata)

for plotkey in plotkeys:
  xlist = list(plotdata[plotkey]["value"])
  ylist = [y for y in plotdata[plotkey]["value"].values()]
  yerr = [y for y in plotdata[plotkey]['variance'].values()]

  plt.errorbar(xlist, ylist, fmt=colors[plotkey], yerr=yerr, label=labels[plotkey])
  
ax = plt.gca()
#ax.set_yscale('log', basey=10) 
plt.title(u'cuboid muscle, serial run')
plt.xlabel('1D elements total')
plt.ylabel('duration (s)')
plt.legend(loc='best')
plt.grid(which='both')
plt.tight_layout()
plt.savefig(output_path+SCENARIO+'_timing_m.png')

######################
# plot times serial
plt.figure(5, figsize=(10,8))
plotdata = dict()
xdata = Set()
plotkeys = Set()

colors = {
  1 : "yo-",
  2 : "ro-",
  3 : "go-",
  4 : "bo-",
  5 : "mo-",
}
labels = {
  1 : "F=1",
  2 : "F=2",
  3 : "F=3",
  4 : "F=4",
  5 : "F=5",
}

for key in datasets:
  
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  f = dataset[6]
  fe = dataset[7]
  user = dataset[14]
  m = dataset[8]
  
  ode_time = dataset[17]
  parabolic_time = dataset[18]
  fe_time = dataset[19]
  pre_fe_time = dataset[20]
  
  ode_time_v = variances[17]
  parabolic_time_v = variances[18]
  fe_time_v = variances[19]
  pre_fe_time_v = variances[20]
  
  plotkey = f
  if plotkey not in plotdata:
    plotdata[plotkey] = dict()
    plotdata[plotkey]['value'] = collections.OrderedDict()
    plotdata[plotkey]['variance'] = collections.OrderedDict()
    
  xvalue = m
  yvalue = parabolic_time
  yvalue_variance = parabolic_time_v
  
  plotdata[plotkey]['value'][xvalue] = yvalue
  plotdata[plotkey]['variance'][xvalue] = yvalue_variance
  xdata.add(xvalue)
  plotkeys.add(plotkey)
  
xlist = list(xdata)

for plotkey in plotkeys:
  xlist = list(plotdata[plotkey]["value"])
  ylist = [y for y in plotdata[plotkey]["value"].values()]
  yerr = [y for y in plotdata[plotkey]['variance'].values()]

  plt.errorbar(xlist, ylist, fmt=colors[plotkey], yerr=yerr, label=labels[plotkey])
  
ax = plt.gca()
#ax.set_yscale('log', basey=10) 
plt.title(u'duration parabolic (1D) solver, cuboid muscle, serial run.\nF is the number of serial half-sarcomeres per muscle fibres')
plt.xlabel('1D elements total')
plt.ylabel('duration (s)')
plt.legend(loc='best')
plt.grid(which='both')
plt.tight_layout()
plt.savefig(output_path+SCENARIO+'_timing_f.png')

if show_plots:
  plt.show()
#quit()

###############################################################
# output to console
print ""
print "------------- duration -------------------------------------------"
print "{:10}, {:6},  {:6}, {:10}, {:10}, {:10}, {:10}".\
format("key", "F", "#M", "init", "stretch", "init", "main")
for key in datasets:
  
  print "{:10}, {:6}, {:6}, {:10}, {:10}, {:10}, {:10}".\
  format(key, datasets[key]["value"][6], datasets[key]["value"][8], 
  fo.str_format_seconds(datasets[key]["value"][10]), 
  fo.str_format_seconds(datasets[key]["value"][11]), 
  fo.str_format_seconds(datasets[key]["value"][12]), 
  fo.str_format_seconds(datasets[key]["value"][13]))
print ""
print ""
  
print "{:10}, {:6}, {:6}, {:10}, {:10}, {:10}, {:10}".\
format("key", "F", "#M", "ODE", "Parabolic", "FE", "pre FE")
for key in datasets:
  
  print "{:10}, {:6}, {:6}, {:10}, {:10}, {:10}, {:10}".\
  format(key, datasets[key]["value"][6], datasets[key]["value"][8], 
  fo.str_format_seconds(datasets[key]["value"][16]), 
  fo.str_format_seconds(datasets[key]["value"][17]), 
  fo.str_format_seconds(datasets[key]["value"][18]), 
  fo.str_format_seconds(datasets[key]["value"][19]))
