import cv2
from brian2 import*
from matplotlib.image import imread
import keyboard


# Put zero to use Laptop webcam or desktop.
# Using one selects secondary webcam
cam = cv2.VideoCapture(0)

tau_M = 10*ms

#### Refactored but taken from Brian2###########
# Basic LIF implementation where the reversal potential
# is replaced with a timed array of the spike events
eqs = '''
dv/dt = (2.2*timing_im(t, i)- v)/tau_M + 0.16*sqrt(2/tau_M)*xi: 1

'''
#################################################

# Functions
#https://docs.opencv.org/3.4/dd/d43/tutorial_py_video_display.html

# Only run loop if the camera can be read corrextly
read_corr = True
# Continually takes image from webcam for live feed
while read_corr:
    # reads input from camera
    # first field boolean if it is read correct - not used
    read_corr, frame = cam.read()
    print(read_corr)
    # resize the image according to fx and fy
    frame = cv2.resize(frame,None,fx=1, fy=0.5)
    # Write the captured image from the webcam to files
    cv2.imwrite("webcam.png",frame)
    # Implement function to convert to event based image[48]
    frame_n = (1-imread('webcam.png'))[::-1, :,0].T
    # Store number of neurons and time in respective variables
    num_samples,N = frame_n.shape
    # Using timed array synatx to store image
    timing_im = TimedArray(frame_n, dt=1*ms)
    im_neuron = NeuronGroup(N, eqs, threshold='v >1', reset='v = 0', method='euler')
    mon = SpikeMonitor(im_neuron)
    mon1 = StateMonitor(im_neuron,'v',record=0)
    # Run brian2 simulation
    run(num_samples*ms,report ='text')
    # Plotting neuron representation of image
    subplot(211)#[48]
    plot(mon.t/ms, mon.i, '.k', ms=3)
    xlim(0,num_samples)
    ylim(0, N)
    xlabel('Time (ms)')
    ylabel('Neuron index')
    subplot(212)
    plot(mon1.t/ms,mon1.v[0])
    # Dont pause whole simulation because of plot "block =False"
    show(block =False)
    # Show the capture to the screen
    cv2.imshow('Input', frame)
    # Hold it on screem
    pause(0.1)
    # Clear the plot to allow next frame to show
    clf()
    # Have exit key to allow for closure 
    if keyboard.is_pressed('esc'):
        break
      
cam.release()



