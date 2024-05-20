import multiprocessing
import pickle


QUEUE_CHUNK_SIZE = 2047


def send_incremental(queue, obj):
    # Serialize the object
    obj_data = pickle.dumps(obj)

    # Calculate the total size
    total_size = len(obj_data)

    # Send the total size first
    queue.put(total_size)

    # Send data in chunks
    for i in range(0, total_size, QUEUE_CHUNK_SIZE):
        # Extract part of the data
        chunk = obj_data[i:i + QUEUE_CHUNK_SIZE]
        queue.put(chunk)


def receive_incremental(queue):
    # Receive the total size first
    total_size = queue.get()

    # Initialize an empty byte array to collect data
    obj_data = bytearray()

    # Receive data until the total size is reached
    while len(obj_data) < total_size:
        # Get the next chunk from the queue
        chunk = queue.get()
        # Extend the bytearray with the new chunk
        obj_data.extend(chunk)

    # Deserialize the object
    return pickle.loads(obj_data)


def _proc_worker(func, queue, args, kwargs):
    try:
        print(f"Running {func.__name__}")
        result = func(*args, **kwargs)
        print(f"Finished {func.__name__}")
        send_incremental(queue, {'result': result, 'error': None})
        print(f"Sent {func.__name__}")
    except Exception as e:
        print(f"Error in {func.__name__}")
        send_incremental(queue, {'result': None, 'error': e})
        print(f"Sent error for {func.__name__}")

class ReturningProcess(multiprocessing.Process):
    def __init__(self, target, args=[], kwargs={}):
        self.target_name = target.__name__
        self.queue = multiprocessing.Queue()
        super().__init__(target=_proc_worker, args=(target, self.queue, args, kwargs))

    def resolve(self):
        print(f"Resolving {self.target_name}")
        val = receive_incremental(self.queue)
        print(f"Received {self.target_name}")
        self.join()
        print(f"Joined {self.target_name}")

        if val['error'] is not None:
            raise val['error']
        else:
            return val['result']
