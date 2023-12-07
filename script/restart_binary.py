from os.path import exists
from os import rename
import os
import struct

def read_line(file):
    bytes_read = file.read(8+8*5+8+8*3)
    if not bytes_read:
        return None
    if len(bytes_read) != 80:
        return None
    return struct.unpack('QdddddNddd', bytes_read)


def read_binary_log(filename):
    data = {key: [] for key in ['step', 'time', 'temperature', 'energy',  'Ea', 'dE', 'selected_atom_id', ] + [f'vac_pos_{i}' for i in range(3)]}

    with open(filename, 'rb') as file:
        while True:
            line = read_line(file)
            if line is None:
                break
            for key, value in zip(data.keys(), line):
                data[key].append(value)
    return data

def write_line(file, line_data):
    # Ensure that line_data has 10 elements (1Q, 5d, 1N, 3d)
    if len(line_data) != 10:
        raise ValueError("Line data must contain 10 elements")

    # Pack the data into a binary format
    packed_data = struct.pack('QdddddNddd', *line_data)
    file.write(packed_data)

def write_binary_log(filename, data, max_steps):
    with open(filename, 'wb') as file:
        for i in range(len(data['step'])):
            if data['step'][i] > max_steps:
                break
            line_data = (
                data['step'][i],
                data['time'][i],
                data['temperature'][i],
                data['energy'][i],
                data['Ea'][i],
                data['dE'][i],
                data['selected_atom_id'][i],
                data['vac_pos_0'][i],
                data['vac_pos_1'][i],
                data['vac_pos_2'][i],
            )
            write_line(file, line_data)

def read_parameters(filename):
    param = {}
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line_split = line.split(' ', 1)
            if len(line_split) == 2:
                param[line_split[0]] = line_split[1].strip('\n')
    return param

if __name__ == '__main__':
    if exists(f'./kmc_log.bin'):
        rename(f'./kmc_log.bin', f'./kmc_log_backup.bin')
    else:
        print("No log file found")
        exit()

    last_step = None
    last_time = None
    last_temperature = None
    last_energy = None

    data=read_binary_log(f'./kmc_log_backup.bin')

    for step, time, temperature, energy in zip(data['step'][::-1], data['time'][::-1], data['temperature'][::-1], data['energy'][::-1]):
        if exists(f'./{step}.cfg.gz'):
            last_step = step
            last_time = time
            last_temperature = temperature
            last_energy = energy
            break
    print(f"{last_step}, {last_time}, {last_temperature}, {last_energy}")
    write_binary_log(f'./kmc_log.bin', data, last_step)

    if exists(f'./kmc_param.txt'):
        rename(f'./kmc_param.txt', f'./kmc_param_backup.txt')

    old_param = read_parameters(f'./kmc_param_backup.txt')
    with open(f'./kmc_param.txt', 'w') as f4:
        f4.write(f"simulation_method {old_param['simulation_method']}\n")
        f4.write(f"json_coefficients_filename {old_param['json_coefficients_filename']}\n")
        if "time_temperature_filename" in old_param.keys():
            f4.write(f"time_temperature_filename {old_param['time_temperature_filename']}\n")
        f4.write(f"config_filename {last_step}.cfg.gz\n")
        f4.write(f"log_dump_steps {old_param['log_dump_steps']}\n")
        f4.write(f"config_dump_steps {old_param['config_dump_steps']}\n")
        f4.write(f"maximum_steps {old_param['maximum_steps']}\n")
        f4.write(f"thermodynamic_averaging_steps {old_param['thermodynamic_averaging_steps']}\n")
        f4.write(f"temperature {last_temperature}\n")
        f4.write(f"element_set {old_param['element_set']}\n")
        f4.write(f"restart_steps {last_step}\n")
        f4.write(f"restart_energy {last_energy}\n")
        f4.write(f"restart_time {last_time}\n")
        f4.write(f"rate_corrector {old_param['rate_corrector']}\n")
        f4.flush()
    print(f"Done...")