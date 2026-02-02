
from .utils import extract_decorator


class Worker(object):
    def __init__(self, input_data, input2, loc):
        self.input_data = input_data
        self.result = None
        self.input2 = input2
        self.loc= loc


    def map(self):
        ff = extract_decorator(self.input_data, self.input2, self.loc)
