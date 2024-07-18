# BEng-Proj
Dissertation Project for BEng.

Project Title: Modelling the behaviour of the Primary Visual Cortex (PVC)

The project uses spiking neural networks. A more biologically realisitc implementation of nodes in ANN's. The project intended to model a singular visual column from the V1 cortex as accurately as possible. It models most neurons from the below diagrams as expanding the architecture proposed unstability in the design which couldnt be solved at the time. The model includes both inhibitory and excitary neuron behaviour. The learning algorithm used is Spike Timing dependent plasticity (STDP) and uses the Izhikevich model to simulate the indvidual neuron behaviour.

![image](https://github.com/user-attachments/assets/1616e305-34a4-4b65-aaa0-6d35404eb55b)
![image](https://github.com/user-attachments/assets/ae7089ab-f2cd-4d28-9af2-2ad89905cee6)

Design structures undsertsood from [1][2]

Included demo of input from webcam, converted into an image represented by neuron spikes to emulate a neuromorphic camera. These inputs then feed into a singular Izhikevich neuron. It is a singular neuron to ensure demonstration ran at acceptable speed.

Example output from model from live camera feed. These emulations of neuromorphic cameras were apart of an exploration of engineering applications of the model, how and where complex asynchronous models could be applied to vision tasks that required power efficiency.
![image](https://github.com/user-attachments/assets/1bb10ed6-982d-42e6-8531-d43bfda69c40)



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
