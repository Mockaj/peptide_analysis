import argparse
import numpy as np

# Set up command line argument parser
parser = argparse.ArgumentParser(description='This scripts filters the sweeps from precycle_sc run with demanded properties, '
'such as specified distance or the angle between the complex peptide helices.')
parser.add_argument('--inputf', '-i', type=str, help='Name of the text file to load. Typically a number of movie files merged into one file.')
parser.add_argument('--outputf', '-o', type=str, help='Name of the text file to write to.')
parser.add_argument('--mind', type=float, help='Minimum distance of peptide center of mass from the membrane center in nanometers.')
parser.add_argument('--maxd', type=float, help='Maximum distance of peptide center of mass from the membrane center in nanometers.')
parser.add_argument('--mina', type=float, help='Minimum angle between peptide helices in degrees.')
parser.add_argument('--maxa', type=float, help='Maximum angle between peptide helices in degrees.')
args = parser.parse_args()

def append_to_file(data):
	with open(args.outputf, 'a') as f:
		for line in data:
			f.write(line)

def calculate_angle(vec1, vec2):
	dot_product = np.dot(vec1, -vec2)
	angle = np.arccos(dot_product)
	# convert radians to degrees
	return np.rad2deg(angle)

def load_input_file(inputf):
		# Load text file
	with open(inputf, 'r') as g:
		print('Preparing for filtering data...')
		lines = g.readlines()
		total_lines = len(lines)
		g.readline()
		g.readline()
		return lines, total_lines

def main():
		lines, total_lines = load_input_file(args.inputf)
		data_holder = ['1538\n', 'sweep 10000; box 17.1351852487 17.1351852487 50.0000000000\n']
		sweep_counter = 0
		mem_cen = 0.0
		pep_cen = 0.0
		molecules_num = 0
		line_counter = 0

		for i in range(0, total_lines):
			line_counter += 1
			percentage = round(line_counter / total_lines * 100)
			if(i%500 == 0):
				print(f'Filtering data {percentage}%. I have accepted {sweep_counter} sweeps so far...', end='\r')
			try:
				molecules_num += 1
				columns = lines[i].split()
				if (molecules_num < 3):
					pep_cen += float(columns[2])
				if (molecules_num == 2):
					pep_cen /= 2
				mem_cen += float(columns[2])
				data_holder.append(lines[i])
				# Angle calculation
				if(int(columns[-1]) == 0):
					if(molecules_num == 1):
						vec1 = np.array([float(columns[3]), float(columns[4]), float(columns[5])])
					elif(molecules_num == 2):
						vec2 = np.array([float(columns[3]), float(columns[4]), float(columns[5])])
						angle = calculate_angle(vec1, vec2)
						angle_is_in_range = True if (angle > args.mina and angle < args.maxa) else False
			# Occurs on the first line of each sweep
			except IndexError:
				mem_cen /= molecules_num
				if(abs(mem_cen-pep_cen) > args.mind and abs(mem_cen-pep_cen) < args.maxd and angle_is_in_range ):
					append_to_file(data_holder)
					sweep_counter += 1
			# Occurs on the second line of each sweep
			except ValueError:
				data_holder = ['1538\n', 'sweep 10000; box 17.1351852487 17.1351852487 50.0000000000\n']
				mem_cen = 0.0
				pep_cen = 0.0
				molecules_num = 0

		print(f'\n number of accepted sweeps = {sweep_counter}')			

if __name__ == "__main__":
	main()