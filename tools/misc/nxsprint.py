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
	#
	# load a TAS nexus file
	#
	def __init__(self, filename):
		file = h5py.File(filename, "r")
		entry = file["entry0"]

		# get scan data
		self.data = entry["data_scan/scanned_variables/data"][:]
		self.data = np.transpose(self.data)
		self.columns = entry["data_scan/scanned_variables/variables_names/property"][:]
		self.columns = np.array([str.decode("utf-8") for str in self.columns])

		self.selected_columns = self.columns

		# get instrument variables
		self.varias = {}
		self.zeros = {}
		self.targets = {}

		instr = entry["instrument"]
		self.instrname = instr["name"][0].decode("utf-8")
		for key in instr.keys():
			varia_path = key + "/value"
			offs_path = key + "/offset_value"
			target_path = key + "/target_value"

			if varia_path in instr:
				self.varias[key] = instr[varia_path][0]
			if offs_path in instr:
				self.zeros[key] = instr[offs_path][0]
			if target_path in instr:
				self.targets[key] = instr[target_path][0]

		# get user infos
		user = entry["user"]
		self.username = user["name"][0].decode("utf-8")
		self.localname = user["namelocalcontact"][0].decode("utf-8")
		self.expnumber = user["proposal"][0].decode("utf-8")

		# get experiment infos
		self.exptitle = entry["title"][0].decode("utf-8")
		self.starttime = entry["start_time"][0].decode("utf-8")
		self.numor = entry["run_number"][0]

		# get sample infos
		sample = entry["sample"]
		self.posqe = (
			sample["qh"][0],
			sample["qk"][0],
			sample["ql"][0],
			sample["en"][0] )
		self.lattice = (
			sample["unit_cell_a"][0],
			sample["unit_cell_b"][0],
			sample["unit_cell_c"][0] )
		self.angles = (
			sample["unit_cell_alpha"][0],
			sample["unit_cell_beta"][0],
			sample["unit_cell_gamma"][0] )
		self.plane0 = ( sample["ax"][0], sample["ay"][0], sample["az"][0] )
		self.plane1 = ( sample["bx"][0], sample["by"][0], sample["bz"][0] )


	#
	# prints a table of the scanned variables
	#
	def print_table(self, table_format = "rounded_grid"):
		indices = np.array([np.where(self.columns == selected_column)[0][0] \
			for selected_column in self.selected_columns])
		print(tab.tabulate(self.data[:,indices], self.columns[indices],
			numalign = "right", tablefmt = table_format))


	#
	# prints the old-style TAS text file format
	#
	def print_retro(self):
		# print header variable
		def print_var(var, name):
			ctr = 0
			for key in self.zeros:
				if ctr % 4 == 0:
					print("\n%s: " % name, end="")
				val = float(self.zeros[key])
				print("%-8s = %6.2f, " % (key, val), end="")
				ctr += 1

		# print header
		print("INSTR: %s" % self.instrname)
		print("USER_: %s" % self.username)
		print("LOCAL: %s" % self.localname)
		print("EXPNO: %s" % self.expnumber)
		print("TITLE: %s" % self.exptitle)
		print("DATE_: %s" % self.starttime)
		print("FILE_: %d" % self.numor)
		print("POSQE: QH = %.4f, QK = %.4f, QL = %.4f, EN = %.4f, UN=meV" % self.posqe)
		print("PARAM: AS = %.4f, BS = %.4f, CS = %.4f" % self.lattice)
		print("PARAM: AA = %.4f, BB = %.4f, CC = %.4f" % self.angles)
		print("PARAM: AX = %.4f, AY = %.4f, AZ = %.4f" % self.plane0)
		print("PARAM: BX = %.4f, BY = %.4f, BZ = %.4f" % self.plane1)
		print_var(self.varias, "VARIA")
		print_var(self.zeros, "ZEROS")
		print_var(self.targets, "TARGET")
		print()

		# print data
		print("DATA_:")
		self.print_table(table_format = "plain")


#
# loads TAS files from the command line and converts them
#
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
