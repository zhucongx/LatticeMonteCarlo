from os.path import exists
from os import rename
import os


def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


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
    if exists(f'./kmc_log.txt'):
        rename(f'./kmc_log.txt', f'./kmc_log_backup.txt')
    last_i = None
    last_step = None
    last_time = None
    last_temperature = None
    last_energy = None
    last_vac_position = None

    for i, line in enumerate(reverse_readline(f'./kmc_log_backup.txt')):
        line = line.split()
        if exists(f'./{line[0]}.cfg.gz'):
            last_i = i
            last_step = line[0]
            last_time = line[1]
            last_temperature = line[2]
            last_energy = line[3]
            # Check if solute_disp data is present (line has 13 columns instead of 10)
            last_vac_position = line[7:10]  # vac1, vac2, vac3
            break
    print(f"{last_i}, {last_step}, {last_time}, {last_temperature}, {last_energy}, {last_vac_position}")
    with open(f'./kmc_log_backup.txt', 'r') as f1, open(f'./kmc_log.txt', 'w') as f2:
        for i, line in enumerate(f1):
            f2.write(line)
            step = line.split()[0]
            if step == last_step:
                f2.flush()
                break
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
        f4.write(f"vacancy_trajectory {last_vac_position[0]} {last_vac_position[1]} {last_vac_position[2]}\n")
        f4.write(f"early_stop {old_param['early_stop']}\n")
        # Add solute_disp parameter if it exists in old parameters
        if "solute_disp" in old_param.keys():
            f4.write(f"solute_disp {old_param['solute_disp']}\n")
        f4.flush()
    print(f"Done...")
