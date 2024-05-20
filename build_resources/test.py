from mp_util import *
from time import sleep
import pickle




def thing():
    print("Hello to thing")
    sleep(2)
    test_obj = 'OMMANIPADMEHUM' * 1000000
    print(len(pickle.dumps(test_obj)))
    return test_obj

foo = ReturningProcess(target=thing)
foo.start()
result = foo.resolve()
print(len(result))
print(result[:10])

