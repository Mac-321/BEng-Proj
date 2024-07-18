# BEng-Proj
Dissertation Project for BEng.

Project Title: Modelling the behaviour of the Primary Visual Cortex (PVC)

The project uses spiking neural networks. A more biologically realisitc implementation of nodes in ANN's. The project intended to model a singular visual column from the V1 cortex as accurately as possible. It models  neurons both inhibitory and excitary neurons. The learning algorithm used is Spike Timing dependent plasticity (STDP) and uses the Izhikevich model to simulate the indvidual neuron behaviour.
![image](https://github.com/user-attachments/assets/bc5e177c-29ee-4887-bf62-0e95d93e4074)![image](https://github.com/user-attachments/assets/0c8d84a6-50e5-4943-ad74-aaeafc1cf753)

Design structures undsertsood from [1][2]

Included demo of input from webcam, converted into an image represented by neuron spikes to emulate a neuromorphic camera. These inputs then feed into a singular Izhikevich neuron. It is a singular neuron to ensure demonstration ran at acceptable speed.


[1] E. R. Kandel, J. H. Schwartz, T. M. Jessell,
S. A. Siegelbaum, and A. J. Hudspeth, Principles
of Neural Science, Fifth Edition. McGraw Hill
Professional, 2012. Accessed: Jan.22, 2023.
[Online]. Available:
https://books.google.co.uk/books/about/Principles
_of_Neural_Science_Fifth_Editi.html?id=Z2yVU
TnlIQsC&redir_esc=y 

[2] S. W. Kuffler, “DISCHARGE PATTERNS
AND FUNCTIONAL ORGANIZATION OF
MAMMALIAN RETINA,” Journal of
Neurophysiology, vol. 16, no. 1, pp. 37–68, Jan.
1953, doi: 10.1152/jn.1953.16.1.37. 
