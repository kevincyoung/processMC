#!/Library/Frameworks/Python.framework/Versions/6.3/Resources/Python.app/Contents/MacOS/Python

from processMC import *

# define a new myField class that defines the number of steps
# and the total length of the control fields

t_final = pi/(4*1.5e6) # seconds
n_steps = 1000
class myField(Field):
	def __init__(self, corrfn = 0.):
		Field.__init__(self, t_final, n_steps, corrfn)

		
# Define Hamiltonian
sig_plus = mat([[0,1],[0,0]])		
sig_minus = sig_plus.H

# Define Hamiltonian
def single_ham(omega_x, omega_z, gamma_x1, gamma_x2, gamma_z1, gamma_z2, epsilon_z, theta):
	h1 = (omega_z + gamma_z1 + gamma_z2 + epsilon_z) * sigZ /2
	h2 = (omega_x + gamma_x1 + gamma_x2 ) * (sig_plus*exp(1.j*theta) + sig_minus*exp(-1.j*theta)) /2
	return h1+h2

gate_name = 'ions_T'
control_field_z = 1.5e6*scipy.ones(n_steps)
control_field_x = scipy.zeros(n_steps)

	
# Add white dephasing noise	
epsilon_z = myField()
epsilon_z.make_white(4.e-4)
epsilon_z.make_noise()

# Add Z control field
omega_z = myField()
omega_z.define_control(control_field_z)
omega_z.make_noise()

# Define white control noise
 # Its amplitude is proportional to the control field
gamma_z1 = myField()
gamma_z1.make_white(4.e-7/((1.5e6)**2))
gamma_z1.define_control(control_field_z)
gamma_z1.make_multiplicative()
gamma_z1.make_noise()

# Define 1/f control noise
 # Its amplitude is proportional to the control field
gamma_z2 = myField()
gamma_z2.make_pink(6.2e2/((1.5e6)**2))
gamma_z2.define_control(control_field_z)
gamma_z2.make_multiplicative()
gamma_z2.make_noise()

# Add X control field
 # Peg control field at maximum value
omega_x = myField()
omega_x.define_control(control_field_x)
omega_x.make_noise()

# Define white control noise
# Its amplitude is proportional to the control field
gamma_x1 = myField()
gamma_x1.make_white(4.e-7/((1.5e6)**2))
gamma_x1.define_control(control_field_x)
gamma_x1.make_multiplicative()
gamma_x1.make_noise()

# Define 1/f control noise
# Its amplitude is proportional to the control field
gamma_x2 = myField()
gamma_x2.make_pink(6.2e2/((1.5e6)**2))
gamma_x2.define_control(control_field_x)
gamma_x2.make_multiplicative()
gamma_x2.make_noise()

# Choose control phase angle
theta = myField()
theta.define_control(scipy.zeros(n_steps))
theta.make_noise()

# Define Liouvillian operator
# ( n_qubits statement isn't strictly necessary )
memphis = Liouvillian(single_ham, omega_x, omega_z, gamma_x1, gamma_x2, gamma_z1, gamma_z2, epsilon_z, theta, verbose = True)
memphis.set_name(gate_name)
tolerance = 1.e-6
memphis.run_converging(tolerance,False,3)
memphis.write_process_matrix(tolerance)

# Propagate the Liouvillian on a parallel machine 10x100 times
# memphis.parallel_propagate(10)

