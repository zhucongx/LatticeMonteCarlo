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


for s in ['s1', 's2', 's3', 's4', 's5']:
    for T in [275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700]:
        if exists(f'{s}/{T}/kmc_log.txt'):
            rename(f'{s}/{T}/kmc_log.txt', f'{s}/{T}/kmc_log_backup.txt')
        last_i = None
        last_step = None
        last_time = None
        last_energy = None

        for i, line in enumerate(reverse_readline(f'{s}/{T}/kmc_log_backup.txt')):
            line = line.split()
            if exists(f'{s}/{T}/map{line[0]}.txt'):
                last_i = i
                last_step = line[0]
                last_time = line[1]
                last_energy = line[2]
                break
        print(f"{s}, {T}K, {last_i}, {last_step}, {last_time}, {last_energy}")
        with open(f'{s}/{T}/kmc_log_backup.txt', 'r') as f1, open(f'{s}/{T}/kmc_log.txt', 'w') as f2:
            for i, line in enumerate(f1):
                f2.write(line)
                step = line.split()[0]
                if step == last_step:
                    f2.flush()
                    break
        with open(f'{s}/{T}/lkmc_param.txt', 'w') as f4:
            f4.write(f"simulation_method ChainKmc\n")
            f4.write(f"map_filename map{last_step}.txt\n")
            f4.write(f"json_coefficients_filename quartic_coefficients.json\n")
            f4.write(f"log_dump_steps 100\n")
            f4.write(f"config_dump_steps 10000\n")
            f4.write(f"maximum_steps 1000000000000000\n")
            f4.write(f"temperature {T}\n")
            f4.write(f"element_set Al Mg Zn\n")
            f4.write(f"restart_steps {last_step}\n")
            f4.write(f"restart_energy {last_energy}\n")
            f4.write(f"restart_time {last_time}\n")
            f4.flush()
        print(f"Done with {s} {T}K")
