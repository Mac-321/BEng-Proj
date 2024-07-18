from brian2 import*
import matplotlib as plt

# Benjamin McCann
# Some code isnt put into functions as it changes the simulation output
# Project uses Brian2 package sourced at [69]

#Converted David Halliday Matlab Function in python [50]
# This function compares individual spike times in the array of input spikes
#to spikes in the output array
def cross_correlation_func(InputSpikes,OutputSpikes,binSize):
    #Largest time difference gives the range of which the
    # cross correlation function will look at
   # number of sampling intervals
    size_array =binSize
    lag_array=[]
    # Assign lag array same number of indexs as the number of intervals
    # e.g. if binsize is 50 it makes an array of length 101 so the graph can
    # be centred at 1 with +/- 50 either side
    for j in range(2*size_array  + 1):
        lag_array.append(0)
    time_lag = 0
    print(len(lag_array))
    neg_tot =50
    lag_tot =101
    neg_offset = neg_tot +0.5
    in_ind =0
    # go through the entire list and then compare each index in input spikes to every value in output spikes
    for i in range(len(InputSpikes)):
        # set back to zero every time loop runs
        out_ind =0
        # Sets end of interval with the Inputspik[i] being -neg_offset from the stop
        # Looks +/- binsize for every index, centred at 1
        #    lag_start < --------- InputSpikes[i] -------> lag_stop
        lag_stop = int(InputSpikes[i]) + neg_offset
        print("current array value",int(InputSpikes[i]))
        lag_start = lag_stop - lag_tot
        print("Lag start value",lag_start)
        print("Lag stop value",lag_stop)
        outputindex =[]
        # Python version of Matlabs "find()" function
        # Will find the index's in the output array
        #that meets the conditions needed.
        # Since the find function finds ALL indexs
        # that meet these requirements in pyhton
        #this is translated to appending an array after going through the array entirely

        # If an index is within this range of -binsize +0.5 -----X----binsize +0.5
        # add to array
        for j in OutputSpikes:
            if(int(j) >lag_start):
                if(int(j) <=lag_stop):
                    outputindex.append(out_ind)   
                    
            out_ind +=1
        print("output index",outputindex)
        
        #go through all values that identify an index in the output
        #array that is a correlated value to the currently viewed input spike time
        # then increment the lag array in the index that correwlates tro the time
        # that the spike pair correlate
        for i in range(len(outputindex)):
                       print("value in output array",int(OutputSpikes[outputindex[i]]))
                       #compare end of range to value in output array
                       indx = math.floor(lag_stop - int(OutputSpikes[outputindex[i]])) + 1
                       print(" Hello" ,indx)
                       # As long as the index is within the array
                       if( indx < len(lag_array)):
                           # Increment the number of synchronous spike pairs
                           # at this time index
                           lag_array[indx] +=1
    plot(lag_array)
    title("Cross-Correlation function")
    xlabel("Lag")
    ylabel("Mean value")
    show()
########### Above from [50]

def plot_output(statm,IVCBETA,IVCALPHA,IVB,L2to3_0,L2to3_1,L2to3_2,L2to3_3,LV_0,LV_1):

    #Plotting input node sim to [46]
    subplot(311)
    plot(statm.t/ms,statm.v[0],label ='Input')
    ylabel('Voltage of spike (mV)')
    legend()
    #Plotting outputs of hidden nodes
    subplot(312)
    plot(LV_1.t/ms,LV_1.v[0],label ='LV_1')
    plot(L2to3_0.t/ms,L2to3_0.v[0],label ='L-2to3[0]')
    plot(L2to3_2.t/ms,L2to3_2.v[0],label ='L-2to3[2]')
    plot(IVCALPHA.t/ms,IVCALPHA.v[0],label ='IVCAlpha')
    plot(IVCBETA.t/ms,IVCBETA.v[0],label ='IVCBeta')
    legend()
    xlabel('Time (ms)')
    ylabel('Voltage of spike (mV)')
    # Plotting outputs of system
    subplot(313)
    plot(IVB.t/ms,IVB.v[0],label ='IVB')
    plot(L2to3_1.t/ms,L2to3_1.v[0],label ='L-2to3[1]')
    plot(L2to3_3.t/ms,L2to3_3.v[0],label ='L-2to3[3]')
    plot(LV_0.t/ms,LV_0.v[0],label ='LV [0]')
    xlabel('Time (ms)')
    ylabel('Voltage of spike (mV)')
    legend()
    tight_layout()
    show()

def plot_synchrony_output(spikmon):
    plot(spikmon.t/ms,spikmon.i,'.')
    # Draw lines at specifc points
    # to help show synchrony
    axvline(x = 50, color = 'r')
    axvline(x = 100, color = 'r')
    axvline(x = 150, color = 'r')
    axvline(x = 200, color = 'r')
    axvline(x = 250, color = 'r')
    axvline(x = 300, color = 'r')
    xlabel("Time (ms)")
    ylabel("Index in PVC model")
    show()

def save_results(statm_1,statm_3,statm_5,statm_7,statm_8):
    ##### Save results!!!!!!!
    # Open or creates a file with the name in the first txt field
    f= open("Results.txt","a+")
    # Write specifics of test
    f.write("\n Parameters e.g.HZ or Amplitude: ")
    f.write("\nInput spike times\n")
    # write the times at which spikes occured in
    #spikemonitor statm_1
    f.write(str(statm_1.t))
    f.write("\nOutput spike times\n ")
    f.write("\n IVB spike time\n")
    f.write(str(statm_3.t))
    f.write("\n L-2to3[1]")
    f.write(str(statm_5.t))
    f.write("\n L-2to3[3]")
    f.write(str(statm_7.t))
    f.write("\n LV [0]") 
    f.write(str(statm_8.t)) 
    f.close()
    
# Neuron membrane relaxation time constant 
tau_M = 10*ms
# Time constants in exponential terms for STDP
#Abbot Song & Miller [45]
tau_pre= 20*ms
tau_post =20*ms
# Node resting potential
v_r = -70
run_time = 700*ms
# Abbot Song & Miller [45] #
gmax = 0.2
dApre = 0.02
dApost = -dApre * tau_pre / tau_post * 1.05
dApost *= gmax
dApre *= gmax
#[45]
stdp_eqn = '''
dApre/dt = -Apre/tau_pre : 1 (event-driven)
dApost/dt = -Apost/tau_post : 1 (event-driven)
w : 1
'''

######
# E.Izhikevich 2006[46]
pre_eqn ='''
I += w
Apre +=  dApre
w = clip(w + Apost, 0, gmax)'''
#############

# Equations used in the Middle/ Output nodes
# Used when using Poisson Input [5][23]
izhikevich_eqns = '''
dv/dt = (0.04*v*v + 5*v + 140 -u + I )/ms + 2*sqrt(2/tau_M)*xi : 1 
du/dt = a*(b*v - u)/ms : 1 
a : 1
b : 1
c : 1
d : 1
I : 1
'''

### Use these equations when using sinusoidal & individual times input only
# Not needed for Poisson input as the syntax works differently
# Sine and individual events stored into an array, this array function of
#time is represented here as I_experiment(t) - then stored as the driving
# current of the node I.
izhikevich_eqns_sine_indv = '''
dv/dt = (0.04*v*v + 5*v + 140 -u + I)/ms +  2*sqrt(2/tau_M)*xi : 1 
du/dt = a*(b*v - u)/ms : 1 
a : 1
b : 1
c : 1
d : 1
I = I_experiment(t): 1
'''

# Use the one below when using Sinusoidal & Individual inputs
input_neuron = NeuronGroup(1,izhikevich_eqns_sine_indv,threshold = 'v >=30',reset ='v = c; u +=d;',method='euler')

# Use the one below when using Poisson inputs
#input_neuron = NeuronGroup(1,izhikevich_eqns,threshold = 'v >=30',reset ='v = c; u +=d;',method='euler')


#------------------------------------------------------------------------
#POISSON INPUT
# Use below for Poisson Input (uncomment line below )
pos = PoissonInput(input_neuron,'v',100,10*Hz,weight=1)
#------------------------------------------------------------------------

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# SINUSOIDAL INPUT
# There are other ways of doing this sinusoidal input, but this syntax in brian2[69]
# was used  because it was needed to input indv events as well
#this timed array syntax was needed
 
# Only run for time 600ms in increments dictated by the system clock
# UNCOMMENT BELOW FOR SINE INPUT
record_t = arange(int(run_time/defaultclock.dt))*defaultclock.dt
# Record the sine wave of Amplitude X with freq Y at times T with spacing
# determined by the default clock (system clock)
# UNCOMMENT BELOW FOR SINE INPUT
I_experiment = TimedArray(15*sin(2*pi*50*Hz*record_t), dt=defaultclock.dt)
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\


# INDIVIDUAL TIMED ARRAY (individual responses)
# Put events into array below - events of different current values
# Stimulates incoming action potential
#UNCOMMENT FOR INDV EVENTS 
#times = [0,6,0,7,0,0,0]
# Events spaced apart by dt value - always starts at zero
#UNCOMMENT FOR INDV EVENTS 
#I_experiment = TimedArray(times,dt = 100*ms)


# Model neuron group
PVC =NeuronGroup(9,izhikevich_eqns,threshold = 'v >=30',reset ='v = c; u +=d;',method='euler')
# All synapses instantiated with unity weight.Trialed random weights (seen in appendix) but caused
# weights to drop to zero and model to stop working.

# Input to IVCBeta 
in_PVC_0 = Synapses(input_neuron,PVC[0],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
in_PVC_0.connect()
in_PVC_0.w =1

# Input to IVCAlpha -correct
in_PVC_1 = Synapses(input_neuron,PVC[1],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
in_PVC_1.connect()
in_PVC_1.w =1


# Input to || & ||| node 1 - correct
in_PVC_4 = Synapses(input_neuron,PVC[4],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
in_PVC_4.connect()
in_PVC_4.w =1


# IVCBeta to IVB - correct
PVC_0_PVC_2 = Synapses(PVC[0],PVC[2],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_0_PVC_2 .connect()
PVC_0_PVC_2 .w =1


# IVCBeta to IVCAlpha
PVC_0_PVC_1 = Synapses(PVC[0],PVC[1],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_0_PVC_1 .connect()
PVC_0_PVC_1 .w =1


# IVCBeta to L || & ||| node 1 - correct 
PVC_0_PVC_4= Synapses(PVC[0],PVC[4],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_0_PVC_4 .connect()
PVC_0_PVC_4.w =1


# IVCAlpha to node 2  L || & ||| - correct
PVC_1_PVC_5= Synapses(PVC[1],PVC[5],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_1_PVC_5 .connect()
PVC_1_PVC_5.w =1


# IVB to node 3 in L || & |||- correct
PVC_2_PVC_6 = Synapses(PVC[2],PVC[6],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_2_PVC_6 .connect()
PVC_2_PVC_6.w =1


# IVB to Node 1 
PVC_2_PVC_5 = Synapses(PVC[2],PVC[5],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_2_PVC_5 .connect()
PVC_2_PVC_5.w =1


# IVB to IVCBeta
PVC_2_PVC_0 = Synapses(PVC[2],PVC[0],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_2_PVC_0 .connect()
PVC_2_PVC_0.w =1

#############################
# Layer || & ||| connections
#############################

# Layer || & ||| node 0 - PVC_3 

# L || & ||| Node 0 to Node 1 - correct
PVC_3_PVC_4 =Synapses(PVC[3],PVC[4],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_3_PVC_4 .connect()
PVC_3_PVC_4.w =1


################
# L || & ||| node 0 to Lv node 1- correxct
PVC_3_PVC_8=Synapses(PVC[3],PVC[8],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_3_PVC_8 .connect()
PVC_3_PVC_8.w =1


# L || & ||| node 0 to Lv node 0 - correct
PVC_3_PVC_7=Synapses(PVC[3],PVC[7],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_3_PVC_7 .connect()
PVC_3_PVC_7.w =1

# node 1 PVC_4
# Node 1 lAYER || & ||| to IVB - correct
PVC4_PVC2 =Synapses(PVC[4],PVC[2],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC4_PVC2 .connect()
PVC4_PVC2.w =1

#
# Node 3 lAYER || & ||| to IVB -corrrect
PVC6_PVC2 =Synapses(PVC[6],PVC[2],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC6_PVC2 .connect()
PVC6_PVC2.w =1

#node 2 - PVC_5
# 2 to 1 - correct
PVC_5_PVC_4 =Synapses(PVC[5],PVC[4],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_5_PVC_4 .connect()
PVC_5_PVC_4.w =1
# Layer || to ||| node 2 to node 0 - correct
PVC_5_PVC_3 =Synapses(PVC[5],PVC[3],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_5_PVC_3 .connect()
PVC_5_PVC_3.w =1

# LV connections - correct
PVC_7_PVC_8  = Synapses(PVC[7],PVC[8],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_7_PVC_8 .connect()
PVC_7_PVC_8.w =1
#
PVC_8_PVC_7  = Synapses(PVC[8],PVC[7],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_8_PVC_7 .connect()
PVC_8_PVC_7.w =1

#
PVC_8_PVC_4  = Synapses(PVC[8],PVC[4],model =stdp_eqn,on_pre=pre_eqn,
             on_post='''Apost += dApost
                     w = clip(w + Apre, 0, gmax)''',method = 'euler')
PVC_8_PVC_4 .connect()
PVC_8_PVC_4.w =1


#Setting Excitatory Regular Firing Patterns for Neurons  [5][23]
#--------------------------------------------------------
input_neuron.a = 0.02
input_neuron.b = 0.2
input_neuron.c = -65
input_neuron.d = 8
input_neuron.u =0.2*-65
input_neuron.v = v_r
#------node 0 L || & |||
PVC[3].a = 0.02
PVC[3].b = 0.2
PVC[3].c = -65
PVC[3].d = 8
PVC[3].u =0.2*-65
PVC[3].v =v_r
#-------------------
PVC[4].a = 0.02
PVC[4].b = 0.2
PVC[4].c = -65
PVC[4].d = 8
PVC[4].u =0.2*-65
PVC[4].v =v_r
#-------------------
PVC[2].a = 0.02
PVC[2].b = 0.2
PVC[2].c = -65
PVC[2].d = 8
PVC[2].u =0.2*-65
PVC[2].v = v_r
#--------------------
PVC[7].a = 0.02
PVC[7].b = 0.2
PVC[7].c = -65
PVC[7].d = 8
PVC[7].u =0.2*-65
PVC[7].v = v_r

#Setting Inhibitatory Thalamo-Cortical Firing Patterns for Neurons [23]
#------------------------------------------------------------------
PVC[1].a = 0.02
PVC[1].b = 0.25
PVC[1].c = -65
PVC[1].d = 0.05
PVC[1].u =0.25*-65
PVC[1].v = v_r
#------------------
PVC[0].a = 0.02
PVC[0].b = 0.25
PVC[0].c = -65
PVC[0].d = 0.05
PVC[0].u =0.25*-65
PVC[0].v = v_r
#-----------------
PVC[5].a = 0.02
PVC[5].b = 0.25
PVC[5].c = -65
PVC[5].d = 0.05
PVC[5].u =0.25*-65
PVC[5].v =v_r
#------------------
PVC[6].a = 0.02
PVC[6].b = 0.25
PVC[6].c = -65
PVC[6].d = 0.05
PVC[6].u =0.25*-65
PVC[6].v =v_r
#-------------------
PVC[8].a = 0.02
PVC[8].b = 0.25
PVC[8].c = -65
PVC[8].d = 0.05
PVC[8].u =0.25*-65
PVC[8].v =v_r
#------------------

#Setting up variables to minitor the states of the different outputs and hidden nodes
statm = StateMonitor(input_neuron,'v', record=0)
statm_1 = SpikeMonitor(input_neuron)
IVCBETA = StateMonitor(PVC[0],'v', record=0)
IVCALPHA = StateMonitor(PVC[1],'v', record=0)
IVB = StateMonitor(PVC[2],'v', record=0)
IVB_SM = SpikeMonitor(PVC[2])
# L|| to ||| --------------------------------
L2to3_0 = StateMonitor(PVC[3],'v', record=0)
L2to3_1 = StateMonitor(PVC[4],'v', record=0)
L2to3_1_SM = SpikeMonitor(PVC[4])
L2to3_2 = StateMonitor(PVC[5],'v', record=0)
L2to3_3 = StateMonitor(PVC[6],'v', record=0)
L2to3_3_SM = SpikeMonitor(PVC[6])
# LV -----------------------------------------
LV_0 = StateMonitor(PVC[7],'v', record=0)
LV_0_SM = SpikeMonitor(PVC[7])
LV_1 = StateMonitor(PVC[8],'v', record=0)

spikmon =SpikeMonitor(PVC)


#Run Experiment
run(run_time,report ='text')

# Plotting Cross-correlation
#cross_correlation_func(statm_1.t/ms,L2to3_1_SM.t/ms,50)

# Save results from test
#save_results(statm_1,IVB_SM,L2to3_1_SM,L2to3_3_SM,LV_0_SM)

# Plot a synchrony plot of spike times of output nodes
#plot_synchrony_output(spikmon)

#Plot input to output response including responses of hidden nodes
plot_output(statm,IVCBETA,IVCALPHA,IVB,L2to3_0,L2to3_1,L2to3_2,L2to3_3,LV_0,LV_1)



