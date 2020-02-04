import pandas as pd

# Read the file(s)
def main():
    x = pd.read_csv('WHUCoV.pack.table', sep='\t', lineterminator='\r')
