import argparse

class ParameterReader:

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-d","--data",required=True,dest="data")

    def readParameters(self):
        arguments = self.parser.parse_args()
        parameters = {}
        parameters["data"] = arguments.data
        return parameters


