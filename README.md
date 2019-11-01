SlippageSim
===========

This tandem repeat simulator was tested with Python 2.7 and BioPython version 1.58.
Further, you need the "hmmemit" program from HMMER version 2 (not 3!)
Get it at http://selab.janelia.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz

To run it try the following commands:

Profile version:
	python2 sim.py <output_directory> <tree_length> <TR_unit_indel_rate_per_unit>

Duplication version:
	python2 simD.py <output_directory> <tree_length> <TR_unit_indel_rate_per_unit>
