import psutil
import multiprocessing as mp
import time


def do_monitoring(outputfile, stop_pipe):
    i = 0
    with open(outputfile, 'w') as out:
        total_mem = psutil.virtual_memory().total
        cpu_count = psutil.cpu_count()
        io_dat = psutil.disk_io_counters()
    
        out.write(f'#{total_mem},{cpu_count}\n')
        while True:
            time.sleep(10)

            cpu_percent = psutil.cpu_percent(interval=3)
            mem_dat = psutil.virtual_memory()

            out.write(f'{time.time()},{cpu_percent},{mem_dat.percent}, {mem_dat.available}, {io_dat.read_time}, {io_dat.write_time}\n')
            if i % 10 == 0:
                out.flush()
            i += 1

            if stop_pipe.poll():
                break



def start_monitor(outputfile):
    stop_getter, stop_sender= mp.Pipe()
    p = mp.Process(target=do_monitoring, args=(outputfile, stop_getter))
    p.start()
    return (p, stop_sender) 

def stop_monitor(stuff):
    p, stop_pipe = stuff
    stop_pipe.send(True)
    p.join(timeout=20)

    if p.is_alive():
        p.terminate()
        p.join(timeout=5)

    print(p.exitcode)

