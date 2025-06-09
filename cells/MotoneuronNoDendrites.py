from .Cell import Cell
from neuron import h

class MotoneuronNoDendrites(Cell):
	"""
	The is a model of the motoneuron soma, as developed by McIntyre 2002.
	This model offers the possibility to simulate the effect of 5-HT as in Booth et al. 1997.
	"""

	def __init__(self,type="WT", drug=True, L=36):
		""" Object initialization.

		Args:
			drug: A boolean flag that is used to decide whether 5-HT is
				inserted in the model or not (default = True).
			L: motoneuron diameter.
		"""
		Cell.__init__(self)

		# Define parameters
		self._drug = drug
		self._L = L
		self._type=type
		self.synapses = []

		self._create_sections()
		self._define_biophysics()

	def _create_sections(self):
		""" Create the sections of the cell. """
		self.soma = h.Section(name='soma',cell=self) # call to cell=self is required to tell NEURON of this object.

	def _define_biophysics(self):
		""" Assign geometry and membrane properties across the cell. """
		self.soma.nseg = 1
		self.soma.cm = 2
		self.soma.Ra = 200
		self.soma.L = self._L
		self.soma.diam = self._L

		self.soma.insert('motoneuron') # Insert the Neuron motoneuron mechanism developed by McIntyre 2002
		if self._drug: self.soma.gcak_motoneuron *= 0.6 #Add the drug effect as in Booth et al 1997

	def current_soma(self, amplitude, duration, delay):
		"""
		Current clamp to motoneuron soma. TODO: remove? 
		"""
		iclamp=h.IClamp(self.soma(0.5))

		iclamp.delay = delay #ms
		iclamp.dur = duration #ms
		iclamp.amp = amplitude #nA

		return iclamp
