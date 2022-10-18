#
# hdf5 tas file converter
# @author Tobias Weber <tweber@ill.fr>
# @date 18-oct-2022
# @license GPLv2
#

import h5py
import numpy as np
import tabulate as tab


class H5Loader:
	def __init__(self, filename):
		file = h5py.File(filename, "r")

		self.data = file["entry0/data_scan/scanned_variables/data"][:]
		self.data = np.transpose(self.data)
		self.columns = file["entry0/data_scan/scanned_variables/variables_names/property"][:]
		self.columns = np.array([str.decode("utf-8") for str in self.columns])

		self.selected_columns = self.columns


	def print_table(self, table_format = "rounded_grid"):
		indices = np.array([np.where(self.columns == selected_column)[0][0] \
			for selected_column in self.selected_columns])
		#indices = [self.columns.index(selected_column) \
		#	for selected_column in self.selected_columns]

		print(tab.tabulate(self.data[:,indices], self.columns[indices],
			numalign = "right", tablefmt = table_format))


	def print_retro(self):
		self.print_table(table_format = "plain")


def main(argv):
	for filename in argv[1:]:
		try:
			h5 = H5Loader(filename)
			h5.selected_columns = [ "QH", "QK", "QL", "EN" ]
			h5.print_retro()
		except FileNotFoundError as err:
			print(err, file = sys.stderr)


if __name__ == "__main__":
	import sys
	main(sys.argv)
