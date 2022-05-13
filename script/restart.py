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


for T in [300, 400, 500, 600, 700, 800]:
    if (exists(f'{T}K/kmc_log.txt')):
        rename(f'{T}K/kmc_log.txt', f'{T}K/kmc_log_backup.txt')
    last_i = None
    last_step = None
    last_time = None
    last_energy = None

    for i, line in enumerate(reverse_readline(f'{T}K/kmc_log_backup.txt')):
        line = line.split()
        if (exists(f'{T}K/{line[0]}.cfg')):
            last_i = i
            last_step = line[0]
            last_time = line[1]
            last_energy = line[2]
            break
    print(f"{T}K, {last_i}, {last_step}, {last_time}, {last_energy}")
    with open(f'{T}K/kmc_log_backup.txt', 'r') as f1, open(f'{T}K/kmc_log.txt', 'w') as f2:
        for i, line in enumerate(f1):
            f2.write(line)
            step = line.split()[0]
            if step == last_step:
                f2.flush()
                break
    with open(f'{T}K/kmc_param.txt', 'w') as f4:
        f4.write(f"config_filename {last_step}.cfg\n")
        f4.write(f"json_coefficients_filename kmc_parameters_quartic.json\n")
        f4.write(f"log_dump_steps 1\n")
        f4.write(f"config_dump_steps 100000\n")
        f4.write(f"maximum_number 10000000000\n")
        f4.write(f"temperature {T}\n")
        f4.write(f"element_symbols Al Mg Zn\n")
        f4.write(f"restart_steps {last_step}\n")
        f4.write(f"restart_energy {last_energy}\n")
        f4.write(f"restart_time {last_time}\n")
        f4.flush()
    print(f"Done with {T}K")
