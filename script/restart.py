from os.path import exists

for T in [300, 600, 900]:
    with open(f'{T}K/kmc_log.txt', 'r') as f1:
        lines = f1.readlines()
    last_i = None
    last_step = None
    last_time = None
    last_energy = None
    for i, line in enumerate(lines):
        line = line.split()
        if exists(f'{T}K/{line[0]}.cfg'):
            last_i = i
            last_step = line[0]
            last_time = line[1]
            last_energy = line[2]

    print(f"{T}K, {last_step}, {last_time}, {last_energy}")

    with open(f'{T}K/kmc_log_backup.txt', 'w') as f2:
        f2.write(''.join(lines))
        f2.flush()
    with open(f'{T}K/kmc_log.txt', 'w') as f3:
        f3.write(''.join(lines[:last_i + 1]))
        f3.flush()
    with open(f'{T}K/kmc.param', 'w') as f4:
        f4.write(f"config_filename {last_step}.cfg\n")
        f4.write(f"json_parameters_filename kmc_parameters_state.json\n")
        f4.write(f"log_dump_steps 1000\n")
        f4.write(f"config_dump_steps 100000\n")
        f4.write(f"maximum_number 10000000000\n")
        f4.write(f"temperature {T}\n")
        f4.write(f"element_symbols Al Mg Zn\n")
        f4.write(f"restart_steps {last_step}\n")
        f4.write(f"restart_energy {last_energy}\n")
        f4.write(f"restart_time {last_time}\n")
        f4.flush()
