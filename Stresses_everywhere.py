#this is the document that basically does combined loading, then will do mohr circle, and then will check the tresca yield criterion
#normal force is zero, #stress due to bending
import numpy as np

moments=np.array([1,2,3,4]) #list of moments across the length of the wing, 
#check if the units fit with compression and tension
max_shear_in__web=np.array([[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]])

height_of_web=np.array([2,3,4,5])

Ixx=np.array([1,2,3,4]) #m^4

    
normal_stresses=moments*height_of_web/2/Ixx   #array containing 

